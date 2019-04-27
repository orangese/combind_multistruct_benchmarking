#!/bin/bash

if [ -z "$1" ] ; then
   TESTS=$(ls -d */tests)
else
   TESTS=$1/tests
fi

echo Running tests in the following directories: $TESTS

$SCHRODINGER/utilities/py.test $TESTS --pypath=$COMBINDHOME
