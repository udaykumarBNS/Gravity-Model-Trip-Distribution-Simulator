# Gravity Model Trip Distribution Simulator

This C++ program implements three types of gravity models used in transportation planning to estimate Origin-Destination (OD) trip matrices based on production, attraction, and inter-zonal distances.

## 🚀 Features

- **Production Constrained Gravity Model**
- **Attraction Constrained Gravity Model**
- **Doubly Constrained Gravity Model**
- Input from user for:
  - Number of zones
  - Production and/or attraction values
  - Distance matrix between zones
- Computes and displays OD trip distribution matrix
- Handles convergence in doubly constrained model using iterative balancing

## 📘 What is a Gravity Model?

Gravity models are used in transportation engineering to estimate the number of trips between zones. These models assume that:
- Trip interactions increase with the size of the zones (in terms of population, jobs, etc.)
- Trip interactions decrease as the distance between zones increases

## 📂 Project Structure

- `main.cpp` — Core implementation of all three gravity models.
- Accepts user input and prints the computed OD matrix.

## 📈 Gravity Models Implemented

### 1. Production Constrained Gravity Model
- Constraints on origin (production) totals.
- Distributes trips from a given zone based on distances and destination attractions.

### 2. Attraction Constrained Gravity Model
- Constraints on destination (attraction) totals.
- Distributes trips arriving at a given zone based on distances and origins.

### 3. Doubly Constrained Gravity Model
- Constraints on both origins and destinations.
- Iteratively adjusts balancing factors (Ai and Bj) for convergence.

## 🧮 Input Format

- Number of zones (integer)
- Production values (for models 1 & 3)
- Attraction values (for models 2 & 3)
- Distance matrix (n x n)
