import numpy as np

arr = np.array([16, 10, 5, 100])
result = np.where(arr == np.amax(arr))

print(result[0][0])

if (result[0][0] == 3):
    print('Last one is the best')
else:
    print('Index: ', result[0][0], ' is the best')

