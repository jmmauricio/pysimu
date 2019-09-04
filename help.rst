python setup.py bdist_wininst

git add .
git commit -m "new commit"
git push origin master

python setup.py sdist
twine upload dist/* -p <password>

hola