.PHONY: install-*

install-uv:
	curl -LsSf https://astral.sh/uv/install.sh | sh

install-cargo:
	curl https://sh.rustup.rs -sSf | sh
