#include <stdio.h>
#include <stdbool.h>

#include <omp.h>

const int data[] = {14, 32, 36, 34, 95, 36, 57, 18, 19, 16};

static inline bool is_even(int x) { return (x & 1) == 0; }

int main()
{
    int sum = 0;

    size_t length = sizeof(data) / sizeof(data[0]);
#pragma omp parallel for reduction (+:sum)
    for (size_t i = 0; i != length; i++)
    {
        if (is_even(data[i]))
        {
            sum += data[i];
        }
    }

    printf("The sum of the even natural numbers is equal to %d.\n", sum);
    return 0;
}