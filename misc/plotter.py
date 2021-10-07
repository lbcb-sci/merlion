#!/usr/bin/env python
import os, sys, argparse, json, seaborn
from matplotlib import pyplot

seaborn.set()
seaborn.set_style("white")
seaborn.despine()
scpb = seaborn.color_palette("Blues")
scpr = seaborn.color_palette("Reds")
scpg = seaborn.cubehelix_palette(rot=-.4)

class Plotter:
  def __init__(self,path, type):
    self.path = path
    self.type = type

  def DrawStack(self, pile):
    if ((self.type == "chimeric" and not pile["is_chimeric_"])):
      return

    figure, ax = pyplot.subplots(1, 1, figsize = (15, 7.5))

    layers = pile["layers_"]
    for i in range(0, len(layers)):
      ax.plot(layers[i]['first'], i + 1, color = scpb[2], marker = 'o', markersize = 2)
      ax.plot(layers[i]['second'], i + 1, color = scpr[2], marker = 'o', markersize = 2)
      ax.plot([layers[i]['first'], layers[i]['second']], [i + 1, i + 1], color = scpg[2], linestyle = ':')

    ax.set_title(pile["id_"])
    figure.text(0.5, 0.05, "base", ha = "center")
    figure.text(0.05, 0.5, "reads", va = "center", rotation = "vertical")
    pyplot.savefig(str(pile["id_"]) + ".png", format = 'png')
    pyplot.close(figure)

  def Run(self):
    try:
      f = open(self.path)
    except Exception:
      print("[merlion::Plotter::Run] error: unable to open file {}".format(self.path))
      return
    try:
      data = json.load(f)
    except Exception:
      print("[merlion::Plotter::Run] error: file is not in JSON format")
      return

    for stack in data:
      self.DrawStack(data[stack])
    return

if __name__ == "__main__":
  parser = argparse.ArgumentParser(
      description = "Plotter is a tool for drawing stacks",
      formatter_class = argparse.ArgumentDefaultsHelpFormatter)
  parser.add_argument("path",
      help = "input file in JSON format")
  parser.add_argument("-t", "--type", default = "all",
      help = "sequence type selection [all, chimeric]")

  args = parser.parse_args()
  plotter = Plotter(args.path, args.type)
  plotter.Run()
