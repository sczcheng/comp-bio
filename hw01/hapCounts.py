import pprint

def readHaps(myData):
    data = open(myData)

    #hapCounts/sample.tsv
    #fin-2.136445439-136692143.tsv

    first_line = data.readline()
    num_tuples = (len(first_line.split()) - 5) * 2 

    final_list = [[] for _ in range(num_tuples)]


    for line in data:
        cells = line.split()[5:]

        for i, cell in enumerate(cells):
            first, second = cell.split(',')
            final_list[i*2].append(first)
            final_list[2*i+1].append(second)

    tuple_list = [tuple(list) for list in final_list]
    return tuple_list

print(readHaps("hapCounts/sample.tsv"))



def hapCounts(hapDataL):
    hapCounts_Dict = {}

    for tup in hapDataL:
        if tup in hapCounts_Dict:
            hapCounts_Dict[tup] += 1
        else:
            hapCounts_Dict[tup] = 1

    return hapCounts_Dict

#pprint.pprint(hapCounts(tuple_list))