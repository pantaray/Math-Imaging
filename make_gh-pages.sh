#!/bin/sh
# Build documentation in `gh-pages` branch

git checkout gh-pages
rm -rf _images _modules _sources _static _stubs
mv df_tools.py df_tools_py.txt
mv loadimg.py loadimg_py.txt
git checkout master docs *.py
git reset HEAD
cd docs/
make html
cd ..
mv -fv docs/build/html/* ./
rm -rf docs *.py *.pyc
mv df_tools_py.txt df_tools.py 
mv loadimg_py.txt loadimg.py
git add -A
git commit -m "Generated gh-pages for `git log master -1 --pretty=short --abbrev-commit`" && git push origin gh-pages ; git checkout master
