from glob import glob
from schrodinger.structure import StructureReader
from schrodinger.structutils.analyze import generate_smiles
import re

pdbs = ['5HK2', '1ZKY', '2XDL', '1BCU', '6FFH', '6DRX', '4MDD', '2XJJ', '1MQ6',
		'2BOH', '4LM0', '2HB8', '5KDT', '1LPZ', '3GBB', '5I2N', '3LPK', '3UDH',
		'3G2N', '5V57', '3B26', '1GWQ', '1GHY', '4LM3', '5T8E', '3B5R', '2XJX',
		'4MRZ', '1KAV', '1W11', '1KSN', '5NFT', '1JVP', '2P2A', '2QRQ', '3GBA',
		'1C1V', '4WIV', '5A0B', '1PD8', '1C1R', '1D3P', '6DK0', '4DDR', '4XNX',
		'3LPI', '2QMG', '3RSX', '3NXV', '3BUG', '1K1I', '1H1S', '5HCV', '3G2K',
		'6DJZ', '4LDO', '3NU0', '1G3D', '4QB3', '1YC4', '1BZJ', '1Z6E', '1E1X',
		'3BUH', '3L59', '3UUO', '3WFG', '1K1J', '4XP4', '2XJG', '1UOM', '1MQJ',
		'5A8Y', '4MSC', '2HAS', '3NTZ', '3RLJ', '1GHW', '4XPA', '1E3G', '1M5E',
		'1MQH', '2Y02', '3VHV', '3S2V', '3RU1', '2QG2', '3AX8', '2XAB', '1GJ7',
		'1KE5', '1XPC', '3A3Z', '4QL8', '3V49', '1EJN', '5KCJ', '4ZO5', '1D3D',
		'1X7E', '2VKM', '2AMB', '1K1N', '3BRA', '6AWP', '6EL7', '4LYW', '5M2V',
		'5I2K', '3B25', '3K22', '3GYF', '3GHW', '1K21', '5TVN', '3CS4', '2J2U',
		'3NYA', '3CKP', '1O5C', '1NL9', '5A0A', '5G5W', '1K1M', '5MWY', '4P6W',
		'4MRW', '1YIM', '6FFI', '4LM2', '2Y04', '1C1U', '3G2H', '5CJ6', '1K1L',
		'1NWL', '5A0C', '2FDP', '2B1V', '2A3I', '3ZYU', '1H1Q', '1L8G', '3L7B',
		'2FAI', '3NXY', '1XKA', '2AA7', '5L7G', '1JSV', '2YI0', '1R5K', '4LKQ',
		'1PYN', '1NFY', '4LLX', '2P4J', '1SQA', '1GHZ', '5ABW', '1NZ7', '1NFU',
		'1SJ0', '1W0Z', '2HAM', '1OWE', '1F0T', '2BQ7', '1E9H', '4XP1', '1TXI',
		'5NEB', '4CSJ', '3NXO', '5I74', '1GJ6', '1C4U', '4NUC', '6FZ4', '6SBT',
		'1C5Z', '2G9Q', '1QXK', '1G36', '5A8E', '5TPA', '2UWD', '2XDK', '4XP5',
		'1F0U', '5CGD', '1BHX', '4QF9', '1GIH', '1GFY', '2PYI', '1GJ8', '1SYH',
		'1YIN', '1MQG', '4LZS', '1FVV', '4LLK', '1ERR', '2BOK', '6DRZ', '1GJC',
		'1IE8', '3U5J', '1YET', '2XHT', '1NO6', '2HAR', '1SQO', '1C83', '2QRH',
		'5V8Q', '1DI8', '1H1R', '1FIN', '1C84', '3NXT', '1OWH', '2QRM', '2G94', 
		'1G7F', '1D4P', '3ZPQ', '3L7A', '3G2L', '1G30', '1XP1', '3KMY', '3B68',
		'4AMJ', '1G32', '6AWO', '1HFP', '2QG0', '1C5T', '3H06', '1BTY', '4DFF',
		'1E1V', '4QIM', '1FPC', '1XP9', '1MQ5', '3AUQ', '1G7G', '2W3A', '4YMB',
		'2XDX', '2QS4', '1S0Z', '1N0T', '2Y03', '1M5C', '1ONZ', '4HF4', '3GHC',
		'3AZ2', '2YCZ', '1C5N', '5H8N', '1GI7', '2QS2', '5MWP', '1NFX', '2YCW',
		'4MR3', '3L7D', '3NY9', '3L79', '4LLP', '1XP6', '5A8X', '3P5O', '1D6W',
		'1GHV', '1WVJ', '1BJV', '4LLJ', '1D9I', '5T8J', '1KMS', '3G2I', '4MSN',
		'5ACY', '4XP6', '2AA5', '1L2I', '2PBW', '5TP9', '2B1Z', '5A09', '3B65',
		'2HB7', '4LM1', '3B0T', '5EQH', '1DM2', '3P0G', '2AM9', '1M4H', '1GI1',
		'3BUF', '4XP9', '4PF3', '2QS3', '2AA2', '1FVT', '3D4S', '3CLD', '4MSA',
		'1CKP', '4LZR', '1C5Y', '1DLR', '6DS0', '4XYA', '1GJA', '1NFW', '4P6X',
		'1H0V', '5L7I', '1T5Z', '1SQT', '1O3P', '1F0S', '1DLS', '1MQI', '3SYR',
		'6DK1', '4N4W', '4Z93', '1ECV', '4NUE', '4OGJ', '3A40', '1G3E', '3NXR',
		'1NNY', '5CGC', '3B66', '3L7C', '1F0R', '6EL9', '4AMI', '1X7R', '1MQD',
		'6DRY', '5L7E', '5A8Z', '1ONY', '2YKI', '2F34', '4LM4', '5NF5', '1FJS',
		'1C88', '5EQI', '1G5S', '2HVC', '3BFT', '1GIJ', '1OHJ', '1C5S', '2YI7',
		'5H8Q', '2QRP', '2QRG', '1LPG', '1KMV', '1G2L', '6EL6', '4UDD', '3E7C',
		'3A2I', '1C5O', '1M5B', '1EB2', '3ZPR', '1O5A', '2Y00', '4KAK', '1K22',
		'3SYM', '1XQC', '3U5L', '3A78', '3KMX', '5NFP', '2AX9', '4LDL', '3I25',
		'1LPK', '1C86', '3AZ3', '5I71', '2F35', '3NY8', '2QFO', '3AZ1', '4XNU',
		'2AYR', '3G2J', '1C5Q', '5G3J', '3B24', '2QS1', '3V4A', '1P93', '2WKY',
		'1S19', '1F5L', '5C1W', '1GZ8', '3B67', '1GJD', '3H03']

pdbs += ['4IB4', '2AXA', '2VT4', '2RH1', '1FKN', '2YEL', '1AQ1', '4M48', '1BOZ',
		 '3Q77', '1A52', '1EZQ', '1A4W', '5EQG', '1YC1', '4OO9', '3BQD', '3WFF',
         '1BJU', '3UI7', '1C5X', '1BZC', '1A8I', '5HK1', '5I6X', '4JKV', '1DB1',
         '1FTM', '1VSO', '5H8H']


raw_files = 'raw_files/*_lig.mae'
pdb_csv = 'pdb.csv'

with open(pdb_csv, 'w') as fp:
	fp.write('"PDB ID","Ligand SMILES"\n')
	for ligand in glob(raw_files):
		
		pdb_id = re.search('raw_files/(.+)\_lig.mae', ligand).group(1)
		pdb_id = pdb_id.upper()

		if pdb_id not in pdbs: continue
		print(ligand)

		with StructureReader(ligand) as st:
			st = list(st)[0]
		smiles = generate_smiles(st)

		fp.write('"{}","{}"\n'.format(pdb_id, smiles))
print('exit.')