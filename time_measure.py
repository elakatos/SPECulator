"""Function to measure full runtime- and sections of the program"""

# Import timer functions
import time

def elapsed_time(previous_time, section_name, cumulative_time):
    current_time = time.time()
    section_time = current_time - previous_time
    cumulative_time += section_time
    print(f"\n{section_name} = {section_time:.3f}s, Total = {cumulative_time:.3f}s\n") # .3f is used to print 3 decimal places
    return current_time, cumulative_time

if __name__ == "__main__":

    start_time = time.time()
    cumulative_time = 0


    # Actual code for the first section
    print("First section")
    # For example, loading data, processing data, etc.
    start_time, cumulative_time_1 = elapsed_time(start_time, "First_section", cumulative_time)

    # Actual code for the second section
    print("Second section")
    start_time, cumulative_time_2 = elapsed_time(start_time, "Second_section", cumulative_time_1)

    print("Third section")
    # More code to execute here
    start_time, cumulative_time_3 = elapsed_time(start_time, "Third_section", cumulative_time_2)

    print("Lasr section")
    # Final bits of code to execute here
    start_time, cumulative_time_4 = elapsed_time(start_time, "Last_section", cumulative_time_3)

    print(f"first: {cumulative_time_1}, second: {cumulative_time_2}, third: {cumulative_time_3}, last: {cumulative_time_4}")