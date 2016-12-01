import random
n = 4
m = 10
r = 10.0
print m, n

# Print c
for i in range(n):
	print round(random.uniform(-r, r), 1),
print

# Print b
for i in range(m):
	print round(random.uniform(-r, r), 1),
print

# Print A
for i in range(n):
	for j in range(m):
		print round(random.uniform(-r, r), 1),
	print
