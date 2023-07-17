
with open("Hp_results.txt") as f:
    list2 = [row.split(',')[0] for row in f]
list2 = list2[1:]
list2 = [float(ll) for ll in list2]
print(list2)