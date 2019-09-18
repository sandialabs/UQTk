#!/bin/sh

find . -type d -print | grep -v ".git" | grep -v doc | grep -v ".\/dep" |grep -v ".\/test"
