#!/usr/bin/env ipython 

from moviepy.editor import *
from moviepy import editor
from moviepy.video.tools.subtitles import SubtitlesClip
import os

# ---------------------- 8< ----------------------------------------------------
def annotate(clip, txt, txt_color= 'grey20', fontsize=50, font='Xolonium-Bold'):
    """ Writes a text at the bottom of the clip. """
    txtclip = editor.TextClip(txt, fontsize=fontsize, font=font, color=txt_color)
    cvc = editor.CompositeVideoClip([clip, txtclip.set_pos(('center', 'bottom'))])
    return cvc.set_duration(clip.duration)
# ---------------------- 8< ----------------------------------------------------

os.chdir("/Users/berald01/Desktop/asciigenome_demo/")
clip = VideoFileClip("bam-3.mov")

sub1= clip.subclip(1.04, 3.17)
sub2= clip.subclip(8.16, 32.00)

#header = TextClip(txt= "ASCIIGenome!\n-    Genome Browser for Terminals    -\n", 
#    font='Amiri-Bold', fontsize=100, bg_color= 'white', color="grey20").set_duration(1)
final= concatenate_videoclips([sub1, sub2], method= 'compose')
final.write_videofile('bam-3.cut.mp4', fps= 3, codec= 'mpeg4')


clip = VideoFileClip("bam-3.cut.mp4")

subs = [((0, 3.17), 'Load bam file'),
        ((3.17, 7.21), 'Go to region'),
        ((7.21, 12), 'Zoom in'),
        ((12, 15), 'Move forward'),
        ((15, 18), 'Zoom out'),
        ((18, 26), 'Filter reads')
        ]
annotated_clips = [annotate(clip.subclip(from_t, to_t), txt, txt_color= 'blue') for (from_t, to_t), txt in subs]
final_clip = editor.concatenate_videoclips(annotated_clips)
final_clip.write_videofile("bam-3.subs.mp4", fps= 3)

# ---------------------- 8< ----------------------------------------------------

os.chdir("/Users/berald01/Desktop/asciigenome_demo/")
clip = VideoFileClip("bigWig-2.mov")

cat= [clip.subclip(20.00, 22.00),
      clip.subclip(24.00, 29.00),
      clip.subclip(30.00, 35.00),
      clip.subclip(41.00, 45),
      clip.subclip(46.5, 48),
      clip.subclip(51, 53),
      clip.subclip(56, 64),
      clip.subclip(80, 82),
]

final= concatenate_videoclips(cat, method= 'compose')
final.write_videofile('bigWig-3.cut.mp4', fps= 3, codec= 'mpeg4')

