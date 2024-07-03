Examples that showcase the usage of PheSA for shape-based docking (negative receptor image) and PheSA screening (flexible and rigid)

The examples can be run by using gradle:

```
./gradlew run
```

The concepts are explained briefly 

Dependencies: PheSA is part of the OpenChemLib  https://github.com/Actelion/openchemlib

PheSA is implemented as a descriptor in the OpenChemLib and requires first a descriptor generation step using a DescriptorHandler. 
There are two DescriptorHandlers for two different purposes: A) generating a PheSA descriptor from a single, bioactive 3D query conformation
and B) generating a descriptor for a candidate molecule, where first a conformer ensemble is generated by the OpenChemLib. 

For case A, the DescriptorHandlerShapeOneConf has to be used, whereas ```StereoMolecule nativeLigandPose```requires 3D coordinates to be present. 

```
DescriptorHandlerShapeOneConf dhsSC = new DescriptorHandlerShapeOneConf();
PheSAMolecule queryShape = dhs.createDescriptor(nativeLigandPose);
```

whereas for case B, DescriptorHandlerShape is used and ```StereoMolecule candidateMol```does not require the presence of 3D Coordinate.

```
DescriptorHandlerShape dhs = new DescriptorHandlerShape();
PheSAMolecule queryShape = dhs.createDescriptor(candidateMol);
```

