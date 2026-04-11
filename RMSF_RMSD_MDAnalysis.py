######RMSD and RMSF analysis through MD Analysis module####

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import MDAnalysis as mda
from MDAnalysis.analysis import rms, align
from MDAnalysisTests.datafiles import PSF,DCD
import warnings


####Understadn the MDA universe###

print(f"Topology file (PSF) :{PSF}")
print(f"Trajectory file (DCD): {DCD}")

u = mda.Universe(PSF, DCD)

print(f"\nSystem information:")
print(f"  Total atoms:     {len(u.atoms)}")
print(f"  Total residues:  {len(u.residues)}")
print(f"  Total segments:  {len(u.segments)}")
print(f"  Trajectory frames: {len(u.trajectory)}")
print(f"  Timestep:        {u.trajectory.dt:.2f} ps")
print(f"  Total time:      {u.trajectory.totaltime:.2f} ps")


####Understadning the atom selection method####

all_atoms = u.select_atoms('all')
ca_atoms = u.select_atoms('name CA')
backbone = u.select_atoms('backbone')
protein = u.select_atoms('protein')
first_resid = u.select_atoms('resid 1-10')
by_resname = u.select_atoms('resname GLY')
heavy_atoms = u.select_atoms('not name H*')
within_5 = u.select_atoms('name CA and around 5.0 resid 50')


print(f"All atoms:              {len(all_atoms)}")
print(f"Cα atoms:               {len(ca_atoms)}")
print(f"Backbone atoms:         {len(backbone)}")
print(f"Protein atoms:          {len(protein)}")
print(f"First 10 residues:      {len(first_resid)}")
print(f"Glycine atoms:          {len(by_resname)}")
print(f"Heavy atoms:            {len(heavy_atoms)}")
print(f"CA within 5Å of res50:  {len(within_5)}")


###Inspect atom properties####
print(f"Frist 5 Calpha atoms:")
for atom in ca_atoms[:5]:
    print(f'Residue {atom.resid} {atom.resname}'
          f"|atom name: {atom.name}"
          f"| position: {atom.position[0]}"
          f"{atom.position[1]}"
          f"{atom.position[2]} Ang")
    
###Aligning before rmsf and rmsd###

aligner = align.AlignTraj(u,u, select='backbone', in_memory=True)
aligner.run()
print("Aligned")

####Calculate RMSF####

rmsf_analysis = rms.RMSF(ca_atoms)
rmsf_analysis.run()
rmsf_mda = rmsf_analysis.results.rmsf

print(f"RMSF computed for {len(rmsf_mda)} Cα atoms")
print(f"Mean RMSF:  {rmsf_mda.mean():.3f} Å")
print(f"Max  RMSF:  {rmsf_mda.max():.3f} Å "
      f"(residue {ca_atoms.resids[np.argmax(rmsf_mda)]})")
print(f"Min  RMSF:  {rmsf_mda.min():.3f} Å "
      f"(residue {ca_atoms.resids[np.argmin(rmsf_mda)]})")


#####Compute RMSD###

rmsd_analysis = rms.RMSD(ca_atoms, reference= None, ref_frame=0)
rmsd_analysis.run()
rmsd_over_time = rmsd_analysis.results.rmsd[:,2]
time_ps = rmsd_analysis.results.rmsd[:,1]

print(f"Mean RMSD over trajectory: {rmsd_over_time.mean()} Å")
print(f"Max  RMSD over trajectory: {rmsd_over_time.max()} Å")

###RMSF adn RMSD plots###
residue_ids = ca_atoms.resids

fig = plt.figure()
gs = gridspec.GridSpec(2,1, figure=fig)


###RMSF over residue###
ax1 = fig.add_subplot(gs[0,0])
ax1.plot(residue_ids, rmsf_mda, color = 'steelblue',linewidth = 1.5, label = 'MDANalsysi RMSF')
ax1.axhline(y=rmsf_mda.mean(), color='gray',
            linestyle='--', alpha=0.7,
            label=f'Mean = {rmsf_mda.mean()} Å')

# Highlight flexible regions (RMSF > mean + std)
threshold = rmsf_mda.mean() + rmsf_mda.std()
flexible  = rmsf_mda > threshold
ax1.fill_between(residue_ids, 0, rmsf_mda,
                 where=flexible,
                 color='red', alpha=0.3,
                 label=f'Flexible regions (>{threshold:.2f} Å)')

ax1.set_xlabel('Residue number', fontsize=12)
ax1.set_ylabel('RMSF (Å)', fontsize=12)
ax1.set_title('Per-residue RMSF — Cα atoms\n'
              '(flexible regions highlighted in red)', fontsize=12)
ax1.legend(fontsize=9)
ax1.grid(True, alpha=0.3)

###RMSD over time####
ax2 = fig.add_subplot(gs[1,0])
ax2.plot(time_ps, rmsd_over_time, color = 'steelblue', linewidth = 1.5)
ax2.axhline(y=rmsd_over_time.mean(), color='red',
            linestyle='--',
            label=f'Mean = {rmsd_over_time.mean():.2f} Å')
ax2.set_xlabel('Time (ps)', fontsize=12)
ax2.set_ylabel('RMSD vs frame 0 (Å)', fontsize=12)
ax2.set_title('RMSD over trajectory', fontsize=12)
ax2.legend()
ax2.grid(True, alpha=0.3)
plt.suptitle('MDAnalysis Trajectory Analysis\n'
             'RMSF and RMSD of Cα atoms', fontsize=14)
plt.tight_layout()
plt.savefig('results_today/w02d4_mdanalysis_trajectory.png', dpi=150)
plt.show()
print("\nSaved w02d4_mdanalysis_trajectory.png")

##Save RMSF data
rmsf_df = pd.DataFrame({
    'residue_id': residue_ids,
    'resname': ca_atoms.resname,
    'rmsf_mda': rmsf_mda

})
rmsf_df.to_csv('results_today/w02d4_rmsf_results.csv', index=False)
print("Saved RMSF data w02d4_rmsf_results.csv")
