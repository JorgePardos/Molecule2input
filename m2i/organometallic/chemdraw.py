"""What RDKit drops when it reads a ChemDraw file: the labels it could not use.

ChemDraw stores a label in one of two ways. A nickname it recognises (CO, Ph,
Ts) is stored together with the atoms it stands for, and RDKit builds them. A
label it does not recognise -- PhCH2OH, or PiPr2 as some people write it -- is
kept as plain text on a node of no particular type, and RDKit turns that node
into a bare carbon. The text is still there, and it is exact, so it is read
here and expanded instead of being lost.

Both flavours of the format are handled: the binary .cdx (a tree of tagged
objects and properties) and the XML .cdxml. Only what is needed is parsed:
which nodes exist, in which order, what element each is, and what it says.
"""

from __future__ import annotations

import re
import struct
from dataclasses import dataclass, field
from pathlib import Path
from xml.etree import ElementTree

from rdkit import Chem

#: CDX object tags that matter here.
FRAGMENT, NODE, TEXT = 0x8003, 0x8004, 0x8006
#: CDX property tags: node type, element, text.
P_TYPE, P_ELEMENT, P_TEXT = 0x0400, 0x0402, 0x0700
#: Node types (CDX numbering) meaning ChemDraw did not work this one out.
UNINTERPRETED = {0, 6, 7}
CONNECTION_POINT = 12
NAMED_TYPES = {
    "unspecified": 0, "element": 1, "elementlist": 2, "elementlistnickname": 3,
    "nickname": 4, "fragment": 5, "formula": 6, "generic": 7,
    "externalconnectionpoint": 12, "linknode": 13,
}
_SYMBOL = re.compile(r"^[A-Z][a-z]?$")


@dataclass
class Object:
    """A CDX object: a fragment, a node, or anything else with children."""

    kind: int
    element: int = 6
    type: int = 1
    text: str = ""
    children: list[Object] = field(default_factory=list)

    def of_kind(self, kind: int) -> list[Object]:
        return [c for c in self.children if c.kind == kind]


@dataclass
class Atom:
    symbol: str
    text: str
    uninterpreted: bool


def labels(path: Path, mol: Chem.Mol) -> dict[int, str]:
    """Atom index in ``mol`` -> the label ChemDraw could not interpret.

    Empty whenever the file cannot be matched to the molecule with certainty:
    a wrong match would put atoms where the chemist never drew them.
    """
    for molecule in molecules(path):
        if [a.symbol for a in molecule] == [a.GetSymbol() for a in mol.GetAtoms()]:
            return {i: a.text for i, a in enumerate(molecule) if a.uninterpreted and a.text}
    return {}


def molecules(path: Path) -> list[list[Atom]]:
    """Every molecule on the page, as the atoms RDKit would build, in order."""
    path = Path(path)
    try:
        if path.suffix.lower() == ".cdx":
            page = _parse_cdx(path.read_bytes())
        elif path.suffix.lower() == ".cdxml":
            page = _parse_cdxml(path.read_text(encoding="utf-8", errors="replace"))
        else:
            return []
    except Exception:  # noqa: BLE001 - a label is a bonus, never a reason to fail
        return []
    return [_atoms(fragment) for fragment in _top_fragments(page)]


# -- the binary format -----------------------------------------------------------


def _parse_cdx(data: bytes) -> Object:
    root = Object(kind=0)
    stack = [root]
    position = 28  # the header
    while position + 2 <= len(data):
        tag = struct.unpack_from("<H", data, position)[0]
        position += 2
        if tag == 0:  # end of the current object
            if len(stack) > 1:
                stack.pop()
            continue
        if tag & 0x8000:
            position += 4  # the object id, which nothing here needs
            child = Object(kind=tag)
            stack[-1].children.append(child)
            stack.append(child)
            continue
        length = struct.unpack_from("<H", data, position)[0]
        position += 2
        if length == 0xFFFF:
            length = struct.unpack_from("<I", data, position)[0]
            position += 4
        value = data[position:position + length]
        position += length
        current = stack[-1]
        if current.kind == NODE and tag == P_TYPE:
            current.type = struct.unpack("<h", value)[0]
        elif current.kind == NODE and tag == P_ELEMENT:
            current.element = struct.unpack("<h", value)[0]
        elif current.kind == TEXT and tag == P_TEXT and len(stack) > 1:
            stack[-2].text = _text(value)
    return root


def _text(value: bytes) -> str:
    runs = struct.unpack_from("<H", value, 0)[0]
    return value[2 + runs * 10:].decode("latin-1", errors="replace").strip()


# -- the XML format --------------------------------------------------------------


def _parse_cdxml(text: str) -> Object:
    return _xml(ElementTree.fromstring(text))


def _xml(element) -> Object:
    kind = {"fragment": FRAGMENT, "n": NODE, "t": TEXT}.get(element.tag, 0)
    node = Object(kind=kind)
    if kind == NODE:
        node.element = int(element.get("Element", 6))
        name = (element.get("NodeType") or "element").lower().replace(" ", "")
        node.type = NAMED_TYPES.get(name, 1)
    for child in element:
        if child.tag == "t":
            node.text = "".join(piece.text or "" for piece in child.iter("s")).strip()
        else:
            node.children.append(_xml(child))
    return node


# -- flattening ------------------------------------------------------------------


def _top_fragments(page: Object) -> list[Object]:
    """The fragments that are molecules, not nickname expansions inside a node."""
    found = []

    def walk(obj: Object) -> None:
        for child in obj.children:
            if child.kind == FRAGMENT:
                found.append(child)
            elif child.kind != NODE:
                walk(child)

    walk(page)
    return found


def _atoms(fragment: Object) -> list[Atom]:
    """The atoms of one fragment, in the order RDKit builds them."""
    out: list[Atom] = []
    for node in fragment.of_kind(NODE):
        expansions = node.of_kind(FRAGMENT)
        if expansions:
            # A nickname: its own atoms take the place of the node.
            for expansion in expansions:
                out += _atoms(expansion)
            continue
        if node.type == CONNECTION_POINT:
            continue
        symbol = Chem.GetPeriodicTable().GetElementSymbol(node.element or 6)
        out.append(Atom(
            symbol=symbol,
            text=node.text,
            uninterpreted=node.type in UNINTERPRETED and bool(node.text)
            and not (_SYMBOL.match(node.text) and node.text == symbol),
        ))
    return out
