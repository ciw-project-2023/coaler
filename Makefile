.DEFAULT_GOAL := release

SHELL := /bin/bash
test: debug
	@echo "Testing..."
	@./build/Debug/test/Test

release:
	@echo "Building in release mode..."
	@source venv/bin/activate && conan build . -s build_type=Release --build missing

debug:
	@echo "Building in debug mode..."
	@source venv/bin/activate && conan build . -s build_type=Debug --build missing

clean:
	@echo "Cleaning..."
	@rm -r ./build/

install:
	@echo "Installing..."
	@install -m755 ./build/Release/src/coaler /usr/local/bin/coaler

container:
	@echo "Building container..."
	@podman build -t coaler:latest .

