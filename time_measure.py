"""Function to measure full runtime- and sections of the program"""

# Import timer functions
import time
from contextlib import contextmanager

@contextmanager
def timer(label):
    start = time.time()
    try:
        yield
    finally:
        end = time.time()
        print(f'{label}: {end - start} seconds')

# Usage:
with timer('Process data'):
    # Your code here
    pass
