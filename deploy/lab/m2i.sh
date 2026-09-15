#!/bin/sh
# Run m2i as a website from this machine. Run from anywhere:
#
#   deploy/lab/m2i.sh check     is this machine ready?
#   deploy/lab/m2i.sh start     build (first time: 15-25 min) and start
#   deploy/lab/m2i.sh url       the public address
#   deploy/lab/m2i.sh status    what is running, and how much memory it uses
#   deploy/lab/m2i.sh logs      follow the logs (Ctrl+C to stop following)
#   deploy/lab/m2i.sh update    pull the latest code, rebuild, restart
#   deploy/lab/m2i.sh stop      stop everything
#
# With TUNNEL_TOKEN in deploy/lab/.env a named tunnel is used (fixed address on
# your domain); without it, a quick tunnel (random trycloudflare.com address).
set -eu

here=$(cd "$(dirname "$0")" && pwd)
script="$here/$(basename "$0")"
root=$(cd "$here/../.." && pwd)
cd "$here"

if [ -f .env ] && grep -q '^TUNNEL_TOKEN=.' .env; then
    profile=named
    service=tunnel
else
    profile=quick
    service=tunnel-quick
fi

compose() {
    docker compose --profile "$profile" "$@"
}

need() {
    command -v "$1" >/dev/null 2>&1 || { echo "missing: $1 -- $2" >&2; return 1; }
}

cmd_check() {
    ok=0
    need docker "install Docker Engine: https://docs.docker.com/engine/install/" || ok=1
    if command -v docker >/dev/null 2>&1; then
        docker compose version >/dev/null 2>&1 \
            || { echo "missing: docker compose plugin (package docker-compose-plugin)" >&2; ok=1; }
        docker info >/dev/null 2>&1 \
            || { echo "docker is installed but not usable by $(id -un): sudo usermod -aG docker $(id -un), then log out and in" >&2; ok=1; }
    fi
    arch=$(uname -m)
    case "$arch" in
        x86_64|aarch64|arm64) echo "architecture  $arch (supported)" ;;
        *) echo "architecture  $arch: not supported (x86_64 or arm64 needed)" >&2; ok=1 ;;
    esac
    memory=$(awk '/MemTotal/ {printf "%d", $2 / 1048576}' /proc/meminfo)
    if [ "$memory" -lt 4 ]; then
        echo "memory        ${memory} GB: at least 4 GB needed (the photo model alone uses 2.4 GB)" >&2; ok=1
    else
        echo "memory        ${memory} GB"
    fi
    disk=$(df -Pk "$root" | awk 'NR == 2 {printf "%d", $4 / 1048576}')
    if [ "$disk" -lt 10 ]; then
        echo "free disk     ${disk} GB: at least 10 GB needed to build the image" >&2; ok=1
    else
        echo "free disk     ${disk} GB"
    fi
    # The tunnel only goes out: to Cloudflare on port 7844. A firewall that
    # blocks it leaves the site unreachable even though everything runs.
    # Only TCP is tested; the tunnels use http2 over TCP, so blocked UDP is fine.
    if ! command -v bash >/dev/null 2>&1 || ! command -v timeout >/dev/null 2>&1; then
        echo "cloudflare    not tested (needs bash and timeout)"
    elif timeout 5 bash -c 'exec 3<>/dev/tcp/region1.v2.argotunnel.com/7844' 2>/dev/null; then
        echo "cloudflare    port 7844 reachable"
    else
        echo "cloudflare    port 7844 NOT reachable: the network blocks the tunnel; ask IT to allow outbound TCP/UDP 7844 to Cloudflare" >&2; ok=1
    fi
    echo "tunnel        $profile"
    [ "$ok" -eq 0 ] && echo "Ready: deploy/lab/m2i.sh start"
    return "$ok"
}

cmd_start() {
    compose up -d --build
    echo "Starting. The address appears in a minute or two: deploy/lab/m2i.sh url"
}

cmd_url() {
    if [ "$profile" = named ]; then
        echo "Named tunnel: the address is the public hostname you set for it in the Cloudflare dashboard."
        return
    fi
    for _ in $(seq 1 30); do
        address=$(compose logs "$service" 2>/dev/null | grep -o 'https://[a-z0-9-]*\.trycloudflare\.com' | tail -n 1 || true)
        if [ -n "$address" ]; then
            echo "$address"
            return
        fi
        sleep 2
    done
    echo "No address yet. Is it running? deploy/lab/m2i.sh status" >&2
    return 1
}

cmd_status() {
    compose ps
    docker stats --no-stream --format 'table {{.Name}}\t{{.MemUsage}}\t{{.CPUPerc}}' $(compose ps -q) 2>/dev/null || true
}

cmd_logs() {
    compose logs -f --tail 100
}

cmd_update() {
    git -C "$root" pull --ff-only
    compose up -d --build
    echo "Updated. With a quick tunnel the address may have changed: deploy/lab/m2i.sh url"
}

cmd_stop() {
    docker compose --profile quick --profile named down
}

case "${1:-}" in
    check|start|url|status|logs|update|stop) "cmd_$1" ;;
    *) sed -n '2,13p' "$script" | sed 's/^# \{0,1\}//'; exit 2 ;;
esac
