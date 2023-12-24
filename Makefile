.DEFAULT_GOAL := build

build:
	@echo "Building..."
	@conan build . -s build_type=Release
