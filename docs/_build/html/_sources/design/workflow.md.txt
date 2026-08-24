# Workflow

I think right now the workflow sould be something like: 

- Define an enum of regions

- define a TrefftMesh depending on those regions

- define a basis

- define a Mapping between those regions and boundary conditions

- create a problem with that mesh[R], basis, R, mapping[R, bcs]

- That way I cannot impose boundary conditions on nonexisting regions and so on.