from os import write

file = open('../Data/cytoBandIdeo.txt').read().splitlines()

final_cytobads = []
# file= file
for line in file:
    if "_" not in line:
        final_cytobads.append(line)

f = open("filteredCytobands.txt", "w")
f.write("\n".join(final_cytobads))

# print(type(file))


