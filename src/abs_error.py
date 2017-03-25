import csv

results = []

def err(actual, expected):
    return abs(expected-actual)

with open("output_fixed.csv", 'r') as csvfile:
    reader = csv.reader(csvfile)
    next(reader, None)
    for row in reader:
        arr = row[:3]
        for method in row[3:]:
            arr.append(round(err(float(row[2]), float(method)), 2))
        results.append(arr)

with open("abs_fixed.csv", 'w') as csvfile:
    writer = csv.writer(csvfile)
    for row in results:
        writer.writerow(row)

results = []

with open("output_floating.csv", 'r') as csvfile:
    reader = csv.reader(csvfile)
    next(reader, None)
    for row in reader:
        arr = row[:4]
        for method in row[4:]:
            arr.append(round(err(float(row[3]), float(method)),2))
        results.append(arr)

with open("abs_floating.csv", 'w') as csvfile:
    writer = csv.writer(csvfile)
    for row in results:
        writer.writerow(row)
