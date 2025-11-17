#!/usr/bin/env python3
import sys

def main():
    if len(sys.argv) != 3:
        print(f"Usage: {sys.argv[0]} n_neurons mpi_size")
        sys.exit(1)

    n_neurons = int(sys.argv[1])
    mpi_size = int(sys.argv[2])

    # 1つ目のループ（小さい方向に探す）
    tmp = mpi_size
    while tmp > 0:
        n_each = (n_neurons + tmp - 1) // tmp
        last_offset = n_each * (tmp - 1)
        if last_offset < n_neurons:
            break
        tmp -= 1
    print(tmp)

    # 2つ目のループ（大きい方向に探す）
    tmp = mpi_size
    while tmp <= 2 * mpi_size:
        n_each = (n_neurons + tmp - 1) // tmp
        last_offset = n_each * (tmp - 1)
        if last_offset < n_neurons:
            break
        tmp += 1
    print(tmp)


if __name__ == "__main__":
    main()
