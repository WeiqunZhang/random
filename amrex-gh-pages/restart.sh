#!/usr/bin/env bash

git clone --single-branch --branch gh-pages git@github.com:AMReX-Codes/amrex.git gh-pages

git clone --mirror git@github.com:AMReX-Codes/amrex.git

cd amrex.git/

git branch -D gh-pages && git reflog expire --expire=now --all && git gc --prune=now --aggressive

git push

cd ..

git clone git@github.com:AMReX-Codes/amrex.git new

cd new/

git checkout --orphan gh-pages

git rm -rf .

cp -rp ../gh-pages/* .

git add .

git commit -m "restart gh-pages"

git push origin gh-pages

cd ..

