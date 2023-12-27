.DEFAULT_GOAL := release

release:
	@echo "Building in release mode..."
	@conan build . -s build_type=Release --build missing

debug:
	@echo "Building in debug mode..."
	@conan build . -s build_type=Debug --build missing

clean:
	@echo "Cleaning..."
	@rm -r ./build/


install:
	@echo "Installing..."
	@install -m 755 ./build/Release/src/coaler /usr/local/bin/coaler