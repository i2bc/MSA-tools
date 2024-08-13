run spectrumany.py
load structure_output.cif
set transparency, 0.4
hide everything
show cartoon, polymer
show sticks, organic or inorganic
select chains_with_occupancy, byres (all and (q > 0))
show surface, chains_with_occupancy
spectrumany q, white yellow red, minimum=0, maximum=100
delete sele chains_with_occupancy
refresh
chains = cmd.get_chains()
for chain in chains: cmd.do(f'extract chain_{chain}, chain {chain}')
group all_chains, chain_*, action=open
delete structure_output
refresh
reset
