# CCG Rotation Suite
These are various reporting scripts that you have to run when you're running production projects.

## Directions for Sarek and its scripts
Please see the top level [README for Sarek](sarek/docs/README.md)

## Install Directions
Please see [Install Directions](./Install.md)

## Running Tests
### Unit Tests
Before releasing this, it should be tested.
1. The following commands can be used to run unit tests and measure unit test coverage.
```bash
cd [git repo dir]
python -m pytest -n auto --cov="." --cov-context=test --cov-report=term-missing tests
wc -l $(find . -name "*.py" -not -path "./scratch/*")
```
Current coverage is 12.57%.