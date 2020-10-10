fin_data = open("hapCounts/sample.tsv")

#fin-2.136445439-136692143.tsv

tuple_list_first = []
tuple_list_second = []
fin_data.readline()

while True:
    s = fin_data.readline()
    
    if s == "":
        break
    
    L = s.split()

    split_list = L[5].split(',')
    #print(split_list)
    tuple_list_first.append(split_list[0])
    tuple_list_second.append(split_list[1])

print(tuple(tuple_list_first))
print(tuple(tuple_list_second))