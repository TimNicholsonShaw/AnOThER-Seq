import distance, time




n = 1000000
list = []
for i in range(n):
    list.append("abcdefghijk")
print ("list made")
start_time = time.time()

count = 0
for item in list:
    distance.hamming('abcdefghijl', item)
    count +=1
    if count%100000==0: print(count/n)



print("--- %s seconds ---" % (time.time() - start_time))