% Homework 3 - Block Diagram for Detector/Decoder
clc; clear; close all;

% Define the block labels
nodes = {
    'Random Bit Generator'
    'Baseband Modulation (Square & Raised Cosine)'
    'Up-Conversion (Multiplication with Cosine)'
    'FFT for Spectral Analysis'
    'Down-Conversion (Multiplication with Cosine)'
    'Filtering using LPF (L=2, L=10)'
    'Threshold Detector (Bit Recovery)'
    'Bit Error Calculation'
};

% Define connections between blocks (directed edges)
edges = [
    1 2  % Random Bit Generator → Baseband Modulation
    2 3  % Baseband Modulation → Up-Conversion
    3 4  % Up-Conversion → FFT for Spectral Analysis
    3 5  % Up-Conversion → Down-Conversion
    5 6  % Down-Conversion → Filtering (L=2, L=10)
    6 7  % Filtering → Threshold Detector
    7 8  % Threshold Detector → Bit Error Calculation
];

% Create directed graph
G = digraph(edges(:,1), edges(:,2), [], nodes);

% Plot the block diagram
figure;
h = plot(G, 'Layout', 'layered', 'NodeFontSize', 12, 'EdgeColor', 'black', ...
    'NodeLabel', nodes, 'ArrowSize', 15, 'LineWidth', 1.5);

% Adjust block (node) appearance
h.NodeLabelColor = 'black';
h.Marker = 's'; % Square nodes for a block-style diagram
h.NodeFontWeight = 'bold';
h.NodeColor = [0.7, 0.85, 1];  % Light blue blocks (valid RGB triplet)

% Title
title('Block Diagram of the Detector/Decoder', 'FontSize', 14, 'FontWeight', 'bold');
