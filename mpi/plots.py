import matplotlib.pyplot as plt

times = [0.012, 0.045, 0.086, 0.519, 1.705]
times_std = [0.003, 0.003, 0.006, 0.076, 0.378]
array_lengths = [10_000, 50_000, 100_000, 500_000, 1_000_000]
plt.errorbar(array_lengths, times, times_std)
plt.xlabel("Array lenght")
plt.ylabel('Times in miliseconds')
plt.title('Send-Recv back and forth an array of ints')
plt.show()