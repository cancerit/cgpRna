#!/usr/bin/env bash
set -e
pytest --cov-branch --cov-report term --cov-report html --cov=run_cgprna --cov-fail-under=10
set +e

# these should not die:

echo -e "\n#################################"
echo      "# Running pycodestyle (style)   #"
echo      "#################################"
pycodestyle run_cgprna

echo -e "\n#########################################"
echo      "# Running radon (cyclomatic complexity) #"
echo      "#########################################"
radon cc -nc run_cgprna

echo -e "\n#########################################"
echo      "# Running radon (maintainability index) #"
echo      "#########################################"
radon mi -s run_cgprna

echo -e "\n##############################"
echo      "# Running mdl (markdownlint) #"
echo      "##############################"
mdl -r ~MD013 .  # ignore line length rule.

exit 0 # don't die based on assements of code quality
