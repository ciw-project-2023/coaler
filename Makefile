.DEFAULT_GOAL := release

test: debug
	@echo "Testing..."
	@./build/Debug/test/Test

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
