source ~/local/gromacs-2021.1-lassen/bin/GMXRC
# 30% cholesterol in a small POPC bilayer - dir cholesterol-inbilers/
python2 ../resources/insane.py -o chol-bilayer.gro -pbc cubic -box 10,10,15 -salt 0.15 -charge auto -l POPC:7 -l CHOL:3 -sol W -p system.top
# update system.top
gmx grompp -p system.top -f ../resources/martini_v2.x_new-rf-em.mdp -c chol-bilayer.gro -maxwarn 10
gmx mdrun -v -c CG-em.gro -nt 8
# make Solvent Rest - gmx make_ndx -f CG-em.gro
gmx grompp -p system.top -f ../resources/martini_v2.x_new-rf-eq1.mdp -c CG-em.gro -maxwarn 10 -n index.ndx
gmx mdrun -v -c CG-eq1.gro -x traj_comp-eq1.xtc -nt 8
gmx grompp -p system.top -f ../resources/martini_v2.x_new-rf-eq4.mdp -c CG-eq1.gro -maxwarn 10 -n index.ndx
gmx mdrun -v -c CG-eq4.gro -x traj_comp-eq4.xtc -nt 8

# x1 cholesterol in vacume - dir cholesterol-single/
# grap inital structure of one chol form (single.gro) from bilayer simulation.
gmx grompp -p system.top -f ../resources/martini_v2.x_new-rf-em.mdp -c single.gro -maxwarn 10
gmx mdrun -v -c CG-em.gro -nt 1
gmx grompp -p system.top -f ../resources/martini_v2.x_new-rf-eq1v.mdp -c CG-em.gro -maxwarn 10
gmx mdrun -v -c CG-eq1.gro -x traj_comp-eq1.xtc -nt 1
gmx grompp -p system.top -f ../resources/martini_v2.x_new-rf-eq4v.mdp -c CG-eq1.gro -maxwarn 10
gmx mdrun -v -c CG-eq4.gro -x traj_comp-eq4.xtc -nt 1

