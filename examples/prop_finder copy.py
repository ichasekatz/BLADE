import matplotlib.pyplot as plt

data = {
1000: [(0.10,-929189.53),(0.15,-929131.51),(0.05,-928972.60),(0.20,-928873.42),(0.25,-928449.97)],
1500: [(0.20,-971214.97),(0.25,-971160.40),(0.15,-971038.66),(0.30,-970904.33),(0.10,-970579.41),(0.35,-970464.93),(0.40,-969853.80),(0.05,-969724.89),(0.45,-969078.16),(0.50,-968141.99)],
2000: [(0.30,-1019190.35),(0.25,-1019085.15),(0.35,-1019062.55),(0.20,-1018722.74),(0.40,-1018717.21),(0.45,-1018163.98),(0.15,-1018063.92),(0.50,-1017408.14),(0.10,-1017039.30),(0.55,-1016451.42)],
2500: [(0.35,-1071936.72),(0.40,-1071908.58),(0.30,-1071701.51),(0.45,-1071629.15),(0.25,-1071183.65),(0.50,-1071105.06),(0.20,-1070352.85),(0.55,-1070338.43),(0.60,-1069327.13),(0.15,-1069160.12)],
3000: [(0.45,-1128929.53),(0.40,-1128890.62),(0.50,-1128681.72),(0.35,-1128557.01),(0.55,-1128149.73),(0.30,-1127914.27),(0.60,-1127331.01),(0.25,-1126939.20),(0.65,-1126217.60),(0.20,-1125595.47)]
}

# plt.figure(figsize=(8,6))

# for T, points in data.items():
#     x=[p[0] for p in points]
#     g=[p[1] for p in points]
#     plt.scatter(x,g,label=f'{T} K')

# plt.xlabel('Hf Composition')
# plt.ylabel('Gibbs Free Energy (J/mol)')
# plt.title('Composition Stability vs Temperature')
# plt.legend()
# plt.tight_layout()
# plt.savefig("/Users/chasekatz/Desktop/School/Research/BLADE/G_vs_x.png", dpi=300)
# plt.show()

plt.figure(figsize=(8,6))

markers = ['o','s','^','D','P']

for i,(T, points) in enumerate(data.items()):
    points = sorted(points, key=lambda p: p[0])

    # convert HF → TA
    x = [1 - p[0] for p in points]   # <-- change here
    g = [p[1] for p in points]

    plt.scatter(x, g, s=70, marker=markers[i], label=f'{T} K')

plt.xlabel('Ta Composition')
plt.ylabel('Gibbs Free Energy (J/mol)')
plt.title('Composition Stability vs Temperature')
plt.xlim(0, 1)
plt.legend()
plt.tight_layout()
plt.savefig("/Users/chasekatz/Desktop/School/Research/BLADE/G_vs_xTa.png", dpi=300)
plt.show()