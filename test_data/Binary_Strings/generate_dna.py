import os
import random

# if not os.path.exists('binary_strings'):
#     os.makedirs('binary_strings')


for i in range(10**6, 10**7, 10**6):
    with open(f"triple_strings{i}_abc.txt", "w") as f:
        f.write(
            "".join(random.choices(["a","b", "c"], k=i))
        )
