fin_data = open("hapCounts/sample.tsv")

#hapCounts/sample.tsv
#fin-2.136445439-136692143.tsv

first_line = fin_data.readline()
num_tuples = (len(first_line.split()) - 5) * 2 

final_list = [[] for _ in range(num_tuples)]


for line in fin_data:
    cells = line.split()[5:]

    for i, cell in enumerate(cells):
        first, second = cell.split(',')
        final_list[i*2].append(first)
        final_list[2*i+1].append(second)

tuple_list = [tuple(list) for list in final_list]

print(tuple_list)
