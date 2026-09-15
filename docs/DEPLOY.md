# Running m2i as a website from your own machine

m2i runs on a Linux computer you already have — a lab workstation, a spare
PC — and is reached from anywhere through a
[Cloudflare Tunnel](https://developers.cloudflare.com/cloudflare-one/connections/connect-networks/):
the machine connects out to Cloudflare, so no port is opened, no public IP is
needed, and the address is HTTPS. Everything is free, photos included.

The site works while the machine is on. The structures people upload are
processed on it and nowhere else (Cloudflare only relays the connection).

## What the machine needs

| | Needed | Checked by `m2i.sh check` |
|---|---|---|
| Linux, x86_64 or arm64 | yes | yes |
| Memory | 4 GB free for m2i (the photo model alone uses 2.4 GB) | total RAM |
| Disk | 10 GB free to build the image | yes |
| Docker Engine with the compose plugin | yes | yes |
| Internet | outbound only: PyPI and Zenodo while building; Cloudflare on port 7844 while running | port 7844 |

No GPU, database or domain is needed.

**Ask first.** Running a public website from a university network may need
the IT department's permission, and some networks block the port the tunnel
uses (outbound 7844). `m2i.sh check` tests the port.

## 1. Install Docker (once)

On Ubuntu or Debian, following [Docker's instructions](https://docs.docker.com/engine/install/ubuntu/):

```bash
curl -fsSL https://get.docker.com | sudo sh
sudo usermod -aG docker $USER
```

Log out and back in, so that your user can run Docker.

## 2. Get m2i onto the machine

The repository is private, so the machine needs access to it — the simplest
is a [personal access token](https://github.com/settings/tokens) with read
access to the repository, used as the password when git asks:

```bash
git clone https://github.com/JorgePardos/Molecule2input.git
cd Molecule2input
```

## 3. Check and start

```bash
sh deploy/lab/m2i.sh check
sh deploy/lab/m2i.sh start
```

The first start builds the image: 15–25 minutes, about 4 GB downloaded
(TensorFlow and the photo model among them). Later starts take seconds.

```bash
sh deploy/lab/m2i.sh url
```

prints the address, something like `https://words-in-a-row.trycloudflare.com`.
Open it, and share it.

The site restarts on its own after a crash or a reboot, as long as Docker
starts at boot (`sudo systemctl enable docker`, the default on most systems).

## Day to day

| | |
|---|---|
| `sh deploy/lab/m2i.sh url` | the current address |
| `sh deploy/lab/m2i.sh status` | what is running and the memory it uses |
| `sh deploy/lab/m2i.sh logs` | follow the logs |
| `sh deploy/lab/m2i.sh update` | pull the latest m2i, rebuild, restart |
| `sh deploy/lab/m2i.sh stop` | stop the site |

From the machine itself, the site is also at http://localhost:7860.

## A fixed address

The free, account-less tunnel above (a *quick tunnel*) gets a **new random
address every time it restarts** — after a reboot or an update — and
Cloudflare describes it as meant for testing: it has no uptime guarantee and a
limit of 200 simultaneous requests, which a research group will not reach.

For an address that never changes you need a domain in a Cloudflare account
(the account is free; a domain costs around 10 € a year, or your institution
may delegate a subdomain):

1. In the Cloudflare dashboard, *Zero Trust → Networks → Tunnels → Create a
   tunnel*, type *Cloudflared*. Give it a name.
2. Copy the token it shows (the long string after `--token`).
3. Under *Public hostnames*, add e.g. `m2i.yourdomain.org`, service
   `HTTP`, URL `app:7860`.
4. On the machine, put the token in `deploy/lab/.env` (git ignores this file):

   ```
   TUNNEL_TOKEN=eyJhIjoi...
   ```

5. `sh deploy/lab/m2i.sh stop`, then `sh deploy/lab/m2i.sh start`.

With a named tunnel you can also restrict who gets in — for example, only
e-mail addresses of your group — with a free
[Cloudflare Access](https://developers.cloudflare.com/cloudflare-one/policies/access/)
policy on that hostname, without changing m2i.

## How it is put together

- `Dockerfile` — the web application (FastAPI serving the API and the page
  on port 7860), DECIMER in its own environment with its weights, run as an
  unprivileged user in hosted mode (`M2I_HOSTED=1`), one worker process,
  uploads capped at 20 MB (`M2I_MAX_UPLOAD_MB`). Sessions are signed cookies;
  every response carries a strict Content-Security-Policy, and the page loads
  nothing from other sites.
- `deploy/lab/compose.yaml` — two containers: `app`, and `cloudflared`, which
  reaches the app over Docker's internal network. The app's port is published
  on `127.0.0.1` only, so the lab network cannot reach it directly. Logs are
  rotated (3 × 10 MB).
- Each visitor works in a private temporary folder; folders untouched for a
  day are removed when new sessions start, so the disk does not fill over
  months.
- The photo model starts loading as soon as the server starts, and stays in
  memory; each photo after that takes a couple of seconds. Whoever opens the
  page in the meantime sees a short "Starting m2i" screen, can already use
  SMILES, ChemDraw files and CIFs, and is let in by the page itself when the
  model is ready.

## Troubleshooting

**`m2i.sh url` prints nothing** — the tunnel waits until the app reports
healthy (up to two minutes after a start). Then `sh deploy/lab/m2i.sh logs`:
lines from `tunnel-quick` saying it cannot connect mean the network blocks
port 7844.

**The address stopped working** — with a quick tunnel, it changed: run
`m2i.sh url` again.

**Cloudflare "Error 1033" although `m2i.sh status` shows everything up** — the
tunnel got an address but cannot connect. Look at
`docker compose -f deploy/lab/compose.yaml --profile quick logs tunnel-quick`:
repeated `failed to dial to edge with quic` means the network drops outbound
UDP. The tunnels here already use `--protocol http2` (TCP 7844) for that
reason; if you run an older copy, `m2i.sh update`.

**"That file is larger than 20 MB"** — the upload limit. Raise it with
`M2I_MAX_UPLOAD_MB` in `deploy/lab/compose.yaml` (under `app: environment:`).

**The machine runs out of memory** — `m2i.sh status` shows who uses it.
m2i needs about 3 GB with the photo model loaded; heavy conformer searches or
metal complexes add to it for a few seconds.

**Building fails while downloading the photo model** — Zenodo, which hosts
the weights, was unreachable. Run `m2i.sh start` again; finished steps are
not repeated.
