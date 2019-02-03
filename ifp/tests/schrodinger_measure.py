from schrodinger.structutils.measure import measure_distance, measure_bond_angle


print("straight line")
print(measure_bond_angle([0, 0, 0], [1, 0, 0], [2, 0, 0]))

print("straight line out of order")
print(measure_bond_angle([0, 0, 0], [2, 0, 0], [1, 0, 0]))

print("acute angle")
print(measure_bond_angle([0, 0, 0], [1, 0, 0], [0, .10, 0]))

print("obtuse angle")
print(measure_bond_angle([0, 0, 0], [1, 0, 0], [2, .10, 0]))

