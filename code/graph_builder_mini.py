import matplotlib.pyplot as plt


s = input()

spots = list(map(float, input().split()))
print(spots)

s_1 = input()

real_parts = list(map(float, input().split()))

s_2 = input()

img_parts = list(map(float, input().split()))

n_spots = int(input())
#spots.pop(n_spots)
#real_parts.pop(n_spots)
#img_parts.pop(n_spots)


plt.figure(figsize=(10, 10))
plt.plot(spots, real_parts, color='blue', marker='x', label="real")
plt.plot(spots, img_parts, color='red', marker='o', label="image")

# Настройки осей и легенды
#plt.yscale('log')        # Логарифмическая ось Y
#plt.ylim(-1, 10**(-5) * 9)
#plt.xlim(spots[0], spots[-1])  # Расширяем пределы по X
plt.xlabel("Ось X")      # Метка оси X
plt.ylabel("Ось Y")      # Метка оси Y
plt.legend()             # Легенда
plt.grid(True)           # Сетка

plt.show()
