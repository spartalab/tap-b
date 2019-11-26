from subprocess import call
net_file = "net/ChicagoRegional_net.txt"
trips_file = "net/ChicagoRegional_trips.txt"


# 2 threads algorithm B
with open('results/2_chicago.txt', 'w') as f:
	call(["./bin/tap", net_file, trips_file, "2"], stdout=f)

# 4 threads algorithm B
with open('results/4_chicago.txt', 'w') as f:
	call(["./bin/tap", net_file, trips_file, "4"], stdout=f)

# 8 threads algorithm B
with open('results/8_chicago.txt', 'w') as f:
	call(["./bin/tap", net_file, trips_file, "8"], stdout=f)

# 16 threads algorithm B
with open('results/16_chicago.txt', 'w') as f:
	call(["./bin/tap", net_file, trips_file, "16"], stdout=f)
