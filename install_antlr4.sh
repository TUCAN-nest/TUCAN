#!/usr/bin/env bash
set -euo pipefail

ANTLR_VERSION=4.13.2
ANTLR_JAR="/usr/local/lib/antlr-${ANTLR_VERSION}-complete.jar"

curl -fsSL "https://www.antlr.org/download/antlr-${ANTLR_VERSION}-complete.jar" -o "${ANTLR_JAR}"

# make antlr4 command available
cat >/usr/local/bin/antlr4 <<'EOF'
#!/usr/bin/env bash
exec java -Xmx500M -jar /usr/local/lib/antlr-4.13.2-complete.jar "$@"
EOF
chmod +x /usr/local/bin/antlr4
