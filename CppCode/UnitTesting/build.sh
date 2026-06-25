#!/bin/bash

# Unit tests
# @source https://stackoverflow.com/questions/13547771/g-project-compilation-with-boost-unit-test
g++ -g TEST_PolygonClass.cpp -o TestPolygon.o -lboost_unit_test_framework -lgmp -lmpfr
g++ -g TEST_GeometryClass.cpp -o TestGeometry.o -lboost_unit_test_framework -lgmp -lmpfr

