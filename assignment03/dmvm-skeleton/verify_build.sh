#!/bin/bash
echo "=== Verification Script ==="
echo ""
echo "1. Source code barrier check:"
grep -n "MPI_Barrier" src/dmvm.c
echo ""
echo "2. Executable timestamp:"
ls -lh exe-ICX
echo ""
echo "3. Source file timestamp:"
ls -lh src/dmvm.c
echo ""
echo "4. Check if barrier in binary:"
strings exe-ICX | grep -i barrier || echo "WARNING: Barrier not found in binary!"
echo ""
echo "5. Active defines in source:"
grep "^#define" src/dmvm.c
