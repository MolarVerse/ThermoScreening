# ThermoScreening

## Description

ThermoScreening is a program to calculate the thermochemical properties of given molecule. The aim is to create a framework to allow screening databases of molecules for their thermochemical properties.

## Development Guide

1. Clone the repository
```bash
git clone https://github.com/stkroe/ThermoScreening.git
```
2. Use git flow
```bash
[master] main
[develop] dev
[version tag prefix] v
```
3. Create a feature branch
```bash
git flow feature start <feature_branch>
```
5. Commit your changes to the feature branch
```bash
git add <files>
git commit -m "message"
git flow feature publish <feature_branch>
```

## TODO
- [ ] Add different engine
- [ ] Add conformer generator
- [ ] Add more tests
- [ ] Add a documentation
- [ ] Add screening framework
