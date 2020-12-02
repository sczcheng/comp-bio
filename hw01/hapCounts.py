import pprint
import collections

def readHaps(file_name):
    data = open(file_name)

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

def x():
    return 3

def hapCounts(hapDataL):
    hapCounts_Dict = collections.defaultdict(int)

    for tup in hapDataL:
        hapCounts_Dict[tup] += 1

    return hapCounts_Dict


tuple_list1 = readHaps("hapCounts/sample.tsv")
print(tuple_list1)
pprint.pprint(collections.Counter(tuple_list1))