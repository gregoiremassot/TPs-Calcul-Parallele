Voici le code réalisé lors des TPs des calcul parallele réalisés avec OpenMPI. Pour compiler sous Ubuntu

* Installer openmpi-bin, libopenmpi-dev
* mpicxx assembling.cc conjugate_gradient.cc grid.cc linear_system.cc main.cc output.cc output_vtk.cc partitioning.cc -o TPs
* mpirun -np 4 TPs