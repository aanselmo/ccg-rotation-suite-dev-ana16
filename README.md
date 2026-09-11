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
/usr/local/miniforge3/bin/coverage run \
    -m pytest tests \
    && /usr/local/miniforge3/bin/coverage report --show-missing
wc -l $(find . -name "*.py" -not -path "./scratch/*" -not -path "*tests*")
```
Current coverage is 9.31%.