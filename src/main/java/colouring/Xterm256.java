package colouring;

import exceptions.InvalidColourException;
import java.awt.Color;
import java.util.LinkedHashMap;
import org.apache.commons.lang3.StringUtils;

public class Xterm256 {

  public static final LinkedHashMap<String, Integer> xtermNameToNumber =
      new LinkedHashMap<String, Integer>();
  public static final LinkedHashMap<Integer, String> contrastColour =
      new LinkedHashMap<Integer, String>();

  // static final HashMap<Integer, String> intColourToName= new HashMap<Integer, String>();

  public Xterm256() {

    // See http://jonasjacek.github.io/colours/
    xtermNameToNumber.put("black", 0);
    xtermNameToNumber.put("maroon", 1);
    xtermNameToNumber.put("green", 2);
    xtermNameToNumber.put("olive", 3);
    xtermNameToNumber.put("navy", 4);
    xtermNameToNumber.put("purple", 5);
    xtermNameToNumber.put("teal", 6);
    xtermNameToNumber.put("silver", 7);
    xtermNameToNumber.put("grey", 8);
    xtermNameToNumber.put("red", 9);
    xtermNameToNumber.put("lime", 10);
    xtermNameToNumber.put("yellow", 11);
    xtermNameToNumber.put("blue", 12);
    xtermNameToNumber.put("fuchsia", 13);
    xtermNameToNumber.put("aqua", 14);
    xtermNameToNumber.put(
        "white",
        231); // NB: white is "15" but 15 shows up as grey-ish. So use 231 (grey100) instead
    xtermNameToNumber.put("grey0", 16);
    xtermNameToNumber.put("navyblue", 17);
    xtermNameToNumber.put("darkblue", 18);
    xtermNameToNumber.put("blue3", 19);
    xtermNameToNumber.put("blue3", 20);
    xtermNameToNumber.put("blue1", 21);
    xtermNameToNumber.put("darkgreen", 22);
    xtermNameToNumber.put("deepskyblue4", 23);
    xtermNameToNumber.put("deepskyblue4", 24);
    xtermNameToNumber.put("deepskyblue4", 25);
    xtermNameToNumber.put("dodgerblue3", 26);
    xtermNameToNumber.put("dodgerblue2", 27);
    xtermNameToNumber.put("green4", 28);
    xtermNameToNumber.put("springgreen4", 29);
    xtermNameToNumber.put("turquoise4", 30);
    xtermNameToNumber.put("deepskyblue3", 31);
    xtermNameToNumber.put("deepskyblue3", 32);
    xtermNameToNumber.put("dodgerblue1", 33);
    xtermNameToNumber.put("green3", 34);
    xtermNameToNumber.put("springgreen3", 35);
    xtermNameToNumber.put("darkcyan", 36);
    xtermNameToNumber.put("lightseagreen", 37);
    xtermNameToNumber.put("deepskyblue2", 38);
    xtermNameToNumber.put("deepskyblue1", 39);
    xtermNameToNumber.put("green3", 40);
    xtermNameToNumber.put("springgreen3", 41);
    xtermNameToNumber.put("springgreen2", 42);
    xtermNameToNumber.put("cyan3", 43);
    xtermNameToNumber.put("darkturquoise", 44);
    xtermNameToNumber.put("turquoise2", 45);
    xtermNameToNumber.put("green1", 46);
    xtermNameToNumber.put("springgreen2", 47);
    xtermNameToNumber.put("springgreen1", 48);
    xtermNameToNumber.put("mediumspringgreen", 49);
    xtermNameToNumber.put("cyan2", 50);
    xtermNameToNumber.put("cyan1", 51);
    xtermNameToNumber.put("darkred", 52);
    xtermNameToNumber.put("deeppink4", 53);
    xtermNameToNumber.put("purple4", 54);
    xtermNameToNumber.put("purple4", 55);
    xtermNameToNumber.put("purple3", 56);
    xtermNameToNumber.put("blueviolet", 57);
    xtermNameToNumber.put("orange4", 58);
    xtermNameToNumber.put("grey37", 59);
    xtermNameToNumber.put("mediumpurple4", 60);
    xtermNameToNumber.put("slateblue3", 61);
    xtermNameToNumber.put("slateblue3", 62);
    xtermNameToNumber.put("royalblue1", 63);
    xtermNameToNumber.put("chartreuse4", 64);
    xtermNameToNumber.put("darkseagreen4", 65);
    xtermNameToNumber.put("paleturquoise4", 66);
    xtermNameToNumber.put("steelblue", 67);
    xtermNameToNumber.put("steelblue3", 68);
    xtermNameToNumber.put("cornflowerblue", 69);
    xtermNameToNumber.put("chartreuse3", 70);
    xtermNameToNumber.put("darkseagreen4", 71);
    xtermNameToNumber.put("cadetblue", 72);
    xtermNameToNumber.put("cadetblue", 73);
    xtermNameToNumber.put("skyblue3", 74);
    xtermNameToNumber.put("steelblue1", 75);
    xtermNameToNumber.put("chartreuse3", 76);
    xtermNameToNumber.put("palegreen3", 77);
    xtermNameToNumber.put("seagreen3", 78);
    xtermNameToNumber.put("aquamarine3", 79);
    xtermNameToNumber.put("mediumturquoise", 80);
    xtermNameToNumber.put("steelblue1", 81);
    xtermNameToNumber.put("chartreuse2", 82);
    xtermNameToNumber.put("seagreen2", 83);
    xtermNameToNumber.put("seagreen1", 84);
    xtermNameToNumber.put("seagreen1", 85);
    xtermNameToNumber.put("aquamarine1", 86);
    xtermNameToNumber.put("darkslategray2", 87);
    xtermNameToNumber.put("darkred", 88);
    xtermNameToNumber.put("deeppink4", 89);
    xtermNameToNumber.put("darkmagenta", 90);
    xtermNameToNumber.put("darkmagenta", 91);
    xtermNameToNumber.put("darkviolet", 92);
    xtermNameToNumber.put("purple", 93);
    xtermNameToNumber.put("orange4", 94);
    xtermNameToNumber.put("lightpink4", 95);
    xtermNameToNumber.put("plum4", 96);
    xtermNameToNumber.put("mediumpurple3", 97);
    xtermNameToNumber.put("mediumpurple3", 98);
    xtermNameToNumber.put("slateblue1", 99);
    xtermNameToNumber.put("yellow4", 100);
    xtermNameToNumber.put("wheat4", 101);
    xtermNameToNumber.put("grey53", 102);
    xtermNameToNumber.put("lightslategrey", 103);
    xtermNameToNumber.put("mediumpurple", 104);
    xtermNameToNumber.put("lightslateblue", 105);
    xtermNameToNumber.put("yellow4", 106);
    xtermNameToNumber.put("darkolivegreen3", 107);
    xtermNameToNumber.put("darkseagreen", 108);
    xtermNameToNumber.put("lightskyblue3", 109);
    xtermNameToNumber.put("lightskyblue3", 110);
    xtermNameToNumber.put("skyblue2", 111);
    xtermNameToNumber.put("chartreuse2", 112);
    xtermNameToNumber.put("darkolivegreen3", 113);
    xtermNameToNumber.put("palegreen3", 114);
    xtermNameToNumber.put("darkseagreen3", 115);
    xtermNameToNumber.put("darkslategray3", 116);
    xtermNameToNumber.put("skyblue1", 117);
    xtermNameToNumber.put("chartreuse1", 118);
    xtermNameToNumber.put("lightgreen", 119);
    xtermNameToNumber.put("lightgreen", 120);
    xtermNameToNumber.put("palegreen1", 121);
    xtermNameToNumber.put("aquamarine1", 122);
    xtermNameToNumber.put("darkslategray1", 123);
    xtermNameToNumber.put("red3", 124);
    xtermNameToNumber.put("deeppink4", 125);
    xtermNameToNumber.put("mediumvioletred", 126);
    xtermNameToNumber.put("magenta3", 127);
    xtermNameToNumber.put("darkviolet", 128);
    xtermNameToNumber.put("purple", 129);
    xtermNameToNumber.put("darkorange3", 130);
    xtermNameToNumber.put("indianred", 131);
    xtermNameToNumber.put("hotpink3", 132);
    xtermNameToNumber.put("mediumorchid3", 133);
    xtermNameToNumber.put("mediumorchid", 134);
    xtermNameToNumber.put("mediumpurple2", 135);
    xtermNameToNumber.put("darkgoldenrod", 136);
    xtermNameToNumber.put("lightsalmon3", 137);
    xtermNameToNumber.put("rosybrown", 138);
    xtermNameToNumber.put("grey63", 139);
    xtermNameToNumber.put("mediumpurple2", 140);
    xtermNameToNumber.put("mediumpurple1", 141);
    xtermNameToNumber.put("gold3", 142);
    xtermNameToNumber.put("darkkhaki", 143);
    xtermNameToNumber.put("navajowhite3", 144);
    xtermNameToNumber.put("grey69", 145);
    xtermNameToNumber.put("lightsteelblue3", 146);
    xtermNameToNumber.put("lightsteelblue", 147);
    xtermNameToNumber.put("yellow3", 148);
    xtermNameToNumber.put("darkolivegreen3", 149);
    xtermNameToNumber.put("darkseagreen3", 150);
    xtermNameToNumber.put("darkseagreen2", 151);
    xtermNameToNumber.put("lightcyan3", 152);
    xtermNameToNumber.put("lightskyblue1", 153);
    xtermNameToNumber.put("greenyellow", 154);
    xtermNameToNumber.put("darkolivegreen2", 155);
    xtermNameToNumber.put("palegreen1", 156);
    xtermNameToNumber.put("darkseagreen2", 157);
    xtermNameToNumber.put("darkseagreen1", 158);
    xtermNameToNumber.put("paleturquoise1", 159);
    xtermNameToNumber.put("red3", 160);
    xtermNameToNumber.put("deeppink3", 161);
    xtermNameToNumber.put("deeppink3", 162);
    xtermNameToNumber.put("magenta3", 163);
    xtermNameToNumber.put("magenta3", 164);
    xtermNameToNumber.put("magenta2", 165);
    xtermNameToNumber.put("darkorange3", 166);
    xtermNameToNumber.put("indianred", 167);
    xtermNameToNumber.put("hotpink3", 168);
    xtermNameToNumber.put("hotpink2", 169);
    xtermNameToNumber.put("orchid", 170);
    xtermNameToNumber.put("mediumorchid1", 171);
    xtermNameToNumber.put("orange3", 172);
    xtermNameToNumber.put("lightsalmon3", 173);
    xtermNameToNumber.put("lightpink3", 174);
    xtermNameToNumber.put("pink3", 175);
    xtermNameToNumber.put("plum3", 176);
    xtermNameToNumber.put("violet", 177);
    xtermNameToNumber.put("gold3", 178);
    xtermNameToNumber.put("lightgoldenrod3", 179);
    xtermNameToNumber.put("tan", 180);
    xtermNameToNumber.put("mistyrose3", 181);
    xtermNameToNumber.put("thistle3", 182);
    xtermNameToNumber.put("plum2", 183);
    xtermNameToNumber.put("yellow3", 184);
    xtermNameToNumber.put("khaki3", 185);
    xtermNameToNumber.put("lightgoldenrod2", 186);
    xtermNameToNumber.put("lightyellow3", 187);
    xtermNameToNumber.put("grey84", 188);
    xtermNameToNumber.put("lightsteelblue1", 189);
    xtermNameToNumber.put("yellow2", 190);
    xtermNameToNumber.put("darkolivegreen1", 191);
    xtermNameToNumber.put("darkolivegreen1", 192);
    xtermNameToNumber.put("darkseagreen1", 193);
    xtermNameToNumber.put("honeydew2", 194);
    xtermNameToNumber.put("lightcyan1", 195);
    xtermNameToNumber.put("red1", 196);
    xtermNameToNumber.put("deeppink2", 197);
    xtermNameToNumber.put("deeppink1", 198);
    xtermNameToNumber.put("deeppink1", 199);
    xtermNameToNumber.put("magenta2", 200);
    xtermNameToNumber.put("magenta1", 201);
    xtermNameToNumber.put("orangered1", 202);
    xtermNameToNumber.put("indianred1", 203);
    xtermNameToNumber.put("indianred1", 204);
    xtermNameToNumber.put("hotpink", 205);
    xtermNameToNumber.put("hotpink", 206);
    xtermNameToNumber.put("mediumorchid1", 207);
    xtermNameToNumber.put("darkorange", 208);
    xtermNameToNumber.put("salmon1", 209);
    xtermNameToNumber.put("lightcoral", 210);
    xtermNameToNumber.put("palevioletred1", 211);
    xtermNameToNumber.put("orchid2", 212);
    xtermNameToNumber.put("orchid1", 213);
    xtermNameToNumber.put("orange1", 214);
    xtermNameToNumber.put("sandybrown", 215);
    xtermNameToNumber.put("lightsalmon1", 216);
    xtermNameToNumber.put("lightpink1", 217);
    xtermNameToNumber.put("pink1", 218);
    xtermNameToNumber.put("plum1", 219);
    xtermNameToNumber.put("gold1", 220);
    xtermNameToNumber.put("lightgoldenrod2", 221);
    xtermNameToNumber.put("lightgoldenrod2", 222);
    xtermNameToNumber.put("navajowhite1", 223);
    xtermNameToNumber.put("mistyrose1", 224);
    xtermNameToNumber.put("thistle1", 225);
    xtermNameToNumber.put("yellow1", 226);
    xtermNameToNumber.put("lightgoldenrod1", 227);
    xtermNameToNumber.put("khaki1", 228);
    xtermNameToNumber.put("wheat1", 229);
    xtermNameToNumber.put("cornsilk1", 230);
    xtermNameToNumber.put("grey100", 231);
    xtermNameToNumber.put("grey3", 232);
    xtermNameToNumber.put("grey7", 233);
    xtermNameToNumber.put("grey11", 234);
    xtermNameToNumber.put("grey15", 235);
    xtermNameToNumber.put("grey19", 236);
    xtermNameToNumber.put("grey23", 237);
    xtermNameToNumber.put("grey27", 238);
    xtermNameToNumber.put("grey30", 239);
    xtermNameToNumber.put("grey35", 240);
    xtermNameToNumber.put("grey39", 241);
    xtermNameToNumber.put("grey42", 242);
    xtermNameToNumber.put("grey46", 243);
    xtermNameToNumber.put("grey50", 244);
    xtermNameToNumber.put("grey54", 245);
    xtermNameToNumber.put("grey58", 246);
    xtermNameToNumber.put("grey62", 247);
    xtermNameToNumber.put("grey66", 248);
    xtermNameToNumber.put("grey70", 249);
    xtermNameToNumber.put("grey74", 250);
    xtermNameToNumber.put("grey78", 251);
    xtermNameToNumber.put("grey82", 252);
    xtermNameToNumber.put("grey85", 253);
    xtermNameToNumber.put("grey89", 254);
    xtermNameToNumber.put("grey93", 255);

    // grey85 and grey19 are the same colours used as default for the fore- and background
    // theme metal. It's better to deal with fewer colours.
    contrastColour.put(0, "grey85");
    contrastColour.put(1, "grey85");
    contrastColour.put(2, "grey85");
    contrastColour.put(3, "grey85");
    contrastColour.put(4, "grey85");
    contrastColour.put(5, "grey85");
    contrastColour.put(6, "grey85");
    contrastColour.put(7, "grey19");
    contrastColour.put(8, "grey85");
    contrastColour.put(9, "grey85");
    contrastColour.put(10, "grey19");
    contrastColour.put(11, "grey19");
    contrastColour.put(12, "grey85");
    contrastColour.put(13, "grey85");
    contrastColour.put(14, "grey19");
    contrastColour.put(15, "grey19");
    contrastColour.put(16, "grey85");
    contrastColour.put(17, "grey85");
    contrastColour.put(18, "grey85");
    contrastColour.put(19, "grey85");
    contrastColour.put(20, "grey85");
    contrastColour.put(21, "grey85");
    contrastColour.put(22, "grey85");
    contrastColour.put(23, "grey85");
    contrastColour.put(24, "grey85");
    contrastColour.put(25, "grey85");
    contrastColour.put(26, "grey85");
    contrastColour.put(27, "grey85");
    contrastColour.put(28, "grey85");
    contrastColour.put(29, "grey85");
    contrastColour.put(30, "grey85");
    contrastColour.put(31, "grey85");
    contrastColour.put(32, "grey85");
    contrastColour.put(33, "grey85");
    contrastColour.put(34, "grey85");
    contrastColour.put(35, "grey85");
    contrastColour.put(36, "grey85");
    contrastColour.put(37, "grey85");
    contrastColour.put(38, "grey85");
    contrastColour.put(39, "grey19");
    contrastColour.put(40, "grey19");
    contrastColour.put(41, "grey19");
    contrastColour.put(42, "grey19");
    contrastColour.put(43, "grey19");
    contrastColour.put(44, "grey19");
    contrastColour.put(45, "grey19");
    contrastColour.put(46, "grey19");
    contrastColour.put(47, "grey19");
    contrastColour.put(48, "grey19");
    contrastColour.put(49, "grey19");
    contrastColour.put(50, "grey19");
    contrastColour.put(51, "grey19");
    contrastColour.put(52, "grey85");
    contrastColour.put(53, "grey85");
    contrastColour.put(54, "grey85");
    contrastColour.put(55, "grey85");
    contrastColour.put(56, "grey85");
    contrastColour.put(57, "grey85");
    contrastColour.put(58, "grey85");
    contrastColour.put(59, "grey85");
    contrastColour.put(60, "grey85");
    contrastColour.put(61, "grey19");
    contrastColour.put(62, "grey19");
    contrastColour.put(63, "grey19");
    contrastColour.put(64, "grey85");
    contrastColour.put(65, "grey85");
    contrastColour.put(66, "grey85");
    contrastColour.put(67, "grey19");
    contrastColour.put(68, "grey19");
    contrastColour.put(69, "grey19");
    contrastColour.put(70, "grey85");
    contrastColour.put(71, "grey19");
    contrastColour.put(72, "grey19");
    contrastColour.put(73, "grey19");
    contrastColour.put(74, "grey19");
    contrastColour.put(75, "grey19");
    contrastColour.put(76, "grey19");
    contrastColour.put(77, "grey19");
    contrastColour.put(78, "grey19");
    contrastColour.put(79, "grey19");
    contrastColour.put(80, "grey19");
    contrastColour.put(81, "grey19");
    contrastColour.put(82, "grey19");
    contrastColour.put(83, "grey19");
    contrastColour.put(84, "grey19");
    contrastColour.put(85, "grey19");
    contrastColour.put(86, "grey19");
    contrastColour.put(87, "grey19");
    contrastColour.put(88, "grey85");
    contrastColour.put(89, "grey85");
    contrastColour.put(90, "grey85");
    contrastColour.put(91, "grey85");
    contrastColour.put(92, "grey85");
    contrastColour.put(93, "grey85");
    contrastColour.put(94, "grey85");
    contrastColour.put(95, "grey85");
    contrastColour.put(96, "grey85");
    contrastColour.put(97, "grey19");
    contrastColour.put(98, "grey19");
    contrastColour.put(99, "grey19");
    contrastColour.put(100, "grey85");
    contrastColour.put(101, "grey85");
    contrastColour.put(102, "grey19");
    contrastColour.put(103, "grey19");
    contrastColour.put(104, "grey19");
    contrastColour.put(105, "grey19");
    contrastColour.put(106, "grey85");
    contrastColour.put(107, "grey19");
    contrastColour.put(108, "grey19");
    contrastColour.put(109, "grey19");
    contrastColour.put(110, "grey19");
    contrastColour.put(111, "grey19");
    contrastColour.put(112, "grey19");
    contrastColour.put(113, "grey19");
    contrastColour.put(114, "grey19");
    contrastColour.put(115, "grey19");
    contrastColour.put(116, "grey19");
    contrastColour.put(117, "grey19");
    contrastColour.put(118, "grey19");
    contrastColour.put(119, "grey19");
    contrastColour.put(120, "grey19");
    contrastColour.put(121, "grey19");
    contrastColour.put(122, "grey19");
    contrastColour.put(123, "grey19");
    contrastColour.put(124, "grey85");
    contrastColour.put(125, "grey85");
    contrastColour.put(126, "grey85");
    contrastColour.put(127, "grey85");
    contrastColour.put(128, "grey85");
    contrastColour.put(129, "grey85");
    contrastColour.put(130, "grey85");
    contrastColour.put(131, "grey19");
    contrastColour.put(132, "grey19");
    contrastColour.put(133, "grey19");
    contrastColour.put(134, "grey19");
    contrastColour.put(135, "grey19");
    contrastColour.put(136, "grey85");
    contrastColour.put(137, "grey19");
    contrastColour.put(138, "grey19");
    contrastColour.put(139, "grey19");
    contrastColour.put(140, "grey19");
    contrastColour.put(141, "grey19");
    contrastColour.put(142, "grey85");
    contrastColour.put(143, "grey19");
    contrastColour.put(144, "grey19");
    contrastColour.put(145, "grey19");
    contrastColour.put(146, "grey19");
    contrastColour.put(147, "grey19");
    contrastColour.put(148, "grey19");
    contrastColour.put(149, "grey19");
    contrastColour.put(150, "grey19");
    contrastColour.put(151, "grey19");
    contrastColour.put(152, "grey19");
    contrastColour.put(153, "grey19");
    contrastColour.put(154, "grey19");
    contrastColour.put(155, "grey19");
    contrastColour.put(156, "grey19");
    contrastColour.put(157, "grey19");
    contrastColour.put(158, "grey19");
    contrastColour.put(159, "grey19");
    contrastColour.put(160, "grey85");
    contrastColour.put(161, "grey85");
    contrastColour.put(162, "grey85");
    contrastColour.put(163, "grey85");
    contrastColour.put(164, "grey85");
    contrastColour.put(165, "grey85");
    contrastColour.put(166, "grey85");
    contrastColour.put(167, "grey19");
    contrastColour.put(168, "grey19");
    contrastColour.put(169, "grey19");
    contrastColour.put(170, "grey19");
    contrastColour.put(171, "grey19");
    contrastColour.put(172, "grey85");
    contrastColour.put(173, "grey19");
    contrastColour.put(174, "grey19");
    contrastColour.put(175, "grey19");
    contrastColour.put(176, "grey19");
    contrastColour.put(177, "grey19");
    contrastColour.put(178, "grey19");
    contrastColour.put(179, "grey19");
    contrastColour.put(180, "grey19");
    contrastColour.put(181, "grey19");
    contrastColour.put(182, "grey19");
    contrastColour.put(183, "grey19");
    contrastColour.put(184, "grey19");
    contrastColour.put(185, "grey19");
    contrastColour.put(186, "grey19");
    contrastColour.put(187, "grey19");
    contrastColour.put(188, "grey19");
    contrastColour.put(189, "grey19");
    contrastColour.put(190, "grey19");
    contrastColour.put(191, "grey19");
    contrastColour.put(192, "grey19");
    contrastColour.put(193, "grey19");
    contrastColour.put(194, "grey19");
    contrastColour.put(195, "grey19");
    contrastColour.put(196, "grey85");
    contrastColour.put(197, "grey85");
    contrastColour.put(198, "grey85");
    contrastColour.put(199, "grey85");
    contrastColour.put(200, "grey85");
    contrastColour.put(201, "grey85");
    contrastColour.put(202, "grey85");
    contrastColour.put(203, "grey19");
    contrastColour.put(204, "grey19");
    contrastColour.put(205, "grey19");
    contrastColour.put(206, "grey19");
    contrastColour.put(207, "grey19");
    contrastColour.put(208, "grey85");
    contrastColour.put(209, "grey19");
    contrastColour.put(210, "grey19");
    contrastColour.put(211, "grey19");
    contrastColour.put(212, "grey19");
    contrastColour.put(213, "grey19");
    contrastColour.put(214, "grey85");
    contrastColour.put(215, "grey19");
    contrastColour.put(216, "grey19");
    contrastColour.put(217, "grey19");
    contrastColour.put(218, "grey19");
    contrastColour.put(219, "grey19");
    contrastColour.put(220, "grey19");
    contrastColour.put(221, "grey19");
    contrastColour.put(222, "grey19");
    contrastColour.put(223, "grey19");
    contrastColour.put(224, "grey19");
    contrastColour.put(225, "grey19");
    contrastColour.put(226, "grey19");
    contrastColour.put(227, "grey19");
    contrastColour.put(228, "grey19");
    contrastColour.put(229, "grey19");
    contrastColour.put(230, "grey19");
    contrastColour.put(231, "grey19");
    contrastColour.put(232, "grey85");
    contrastColour.put(233, "grey85");
    contrastColour.put(234, "grey85");
    contrastColour.put(235, "grey85");
    contrastColour.put(236, "grey85");
    contrastColour.put(237, "grey85");
    contrastColour.put(238, "grey85");
    contrastColour.put(239, "grey85");
    contrastColour.put(240, "grey85");
    contrastColour.put(241, "grey85");
    contrastColour.put(242, "grey85");
    contrastColour.put(243, "grey85");
    contrastColour.put(244, "grey19");
    contrastColour.put(245, "grey19");
    contrastColour.put(246, "grey19");
    contrastColour.put(247, "grey19");
    contrastColour.put(248, "grey19");
    contrastColour.put(249, "grey19");
    contrastColour.put(250, "grey19");
    contrastColour.put(251, "grey19");
    contrastColour.put(252, "grey19");
    contrastColour.put(253, "grey19");
    contrastColour.put(254, "grey19");
    contrastColour.put(255, "grey19");
  }

  public static final Color xterm256ToColour(int xterm256) throws InvalidColourException {

    String rgb = xterm256toRGB(xterm256);
    int r = Integer.parseInt(rgb.substring(0, 2), 16);
    int g = Integer.parseInt(rgb.substring(2, 4), 16);
    int b = Integer.parseInt(rgb.substring(4, 6), 16);

    float[] hsb = Color.RGBtoHSB(r, g, b, null);
    return Color.getHSBColor(hsb[0], hsb[1], hsb[2]);
  }

  //	public static final LinkedHashMap<String, Integer> mapColourNameToXterm256(){
  //		return xtermNameToNumber;
  //	}
  //	public static final LinkedHashMap<Integer, String> contrastColourMap(){
  //		return contrastColour;
  //	}
  /**
   * Return the colour integer for the given input string. colourName can be: 1) The name itself
   * like "red" 2) An integer string in the range 0-255, e.g. "230" 3) Prefix of a name, e.g.
   * "lightsal". Colour names are case insensitive.
   */
  public static int colourNameToXterm256(String colourName) throws InvalidColourException {

    try { // See if colour name is already an int between 0 and 255
      int xterm = Integer.parseInt(colourName);
      if (xterm >= 0 && xterm <= 255) {
        return xterm;
      }
    } catch (NumberFormatException e) {
      //
    }

    colourName = colourName.toLowerCase();

    if (xtermNameToNumber.containsKey(colourName)) {
      return xtermNameToNumber.get(colourName);
    } else {
      // Try approximate matching
      for (String x : xtermNameToNumber.keySet()) {
        if (x.startsWith(colourName)) {
          return xtermNameToNumber.get(x);
        }
      }
      System.err.println("Unrecognized colour: " + colourName);
      throw new InvalidColourException();
    }
  }

  public static String colourShowForTerminal() throws InvalidColourException {

    int maxLen = 0;
    for (String x : xtermNameToNumber.keySet()) {
      if (x.length() > maxLen) {
        maxLen = x.length();
      }
    }

    StringBuilder sb = new StringBuilder();
    int i = 0;
    for (String x : xtermNameToNumber.keySet()) {
      i++;
      int xterm = xtermNameToNumber.get(x);
      int spacer = maxLen - x.length();
      sb.append(
          xterm
              + ": \033[38;5;"
              + xterm
              + "m"
              + x
              + "\033[38;5;"
              + Config.get256Colour(ConfigKey.foreground)
              + ";48;5;"
              + Config.get256Colour(ConfigKey.background)
              + "m");

      if (i == 3) { // Arrange colours in this many columns
        sb.append("\n");
        i = 0;
      } else {
        sb.append(StringUtils.repeat(" ", spacer + 1));
      }
    }
    return sb.toString();
  }

  private static final String xterm256toRGB(int xterm256) throws InvalidColourException {

    if (xterm256 < 0 || xterm256 > 255) {
      throw new InvalidColourException();
    }

    // From https://gist.github.com/MicahElliott/719710

    // Primary 3-bit (8 colours). Unique representation!
    if (xterm256 == 0) {
      return "000000";
    }
    if (xterm256 == 1) {
      return "800000";
    }
    if (xterm256 == 2) {
      return "008000";
    }
    if (xterm256 == 3) {
      return "808000";
    }
    if (xterm256 == 4) {
      return "000080";
    }
    if (xterm256 == 5) {
      return "800080";
    }
    if (xterm256 == 6) {
      return "008080";
    }
    if (xterm256 == 7) {
      return "c0c0c0";
    }

    // Equivalent "bright" versions of original 8 colours.
    if (xterm256 == 8) {
      return "808080";
    }
    if (xterm256 == 9) {
      return "ff0000";
    }
    if (xterm256 == 10) {
      return "00ff00";
    }
    if (xterm256 == 11) {
      return "ffff00";
    }
    if (xterm256 == 12) {
      return "0000ff";
    }
    if (xterm256 == 13) {
      return "ff00ff";
    }
    if (xterm256 == 14) {
      return "00ffff";
    }
    if (xterm256 == 15) {
      return "ffffff";
    }

    // Strictly ascending.
    if (xterm256 == 16) {
      return "000000";
    }
    if (xterm256 == 17) {
      return "00005f";
    }
    if (xterm256 == 18) {
      return "000087";
    }
    if (xterm256 == 19) {
      return "0000af";
    }
    if (xterm256 == 20) {
      return "0000d7";
    }
    if (xterm256 == 21) {
      return "0000ff";
    }
    if (xterm256 == 22) {
      return "005f00";
    }
    if (xterm256 == 23) {
      return "005f5f";
    }
    if (xterm256 == 24) {
      return "005f87";
    }
    if (xterm256 == 25) {
      return "005faf";
    }
    if (xterm256 == 26) {
      return "005fd7";
    }
    if (xterm256 == 27) {
      return "005fff";
    }
    if (xterm256 == 28) {
      return "008700";
    }
    if (xterm256 == 29) {
      return "00875f";
    }
    if (xterm256 == 30) {
      return "008787";
    }
    if (xterm256 == 31) {
      return "0087af";
    }
    if (xterm256 == 32) {
      return "0087d7";
    }
    if (xterm256 == 33) {
      return "0087ff";
    }
    if (xterm256 == 34) {
      return "00af00";
    }
    if (xterm256 == 35) {
      return "00af5f";
    }
    if (xterm256 == 36) {
      return "00af87";
    }
    if (xterm256 == 37) {
      return "00afaf";
    }
    if (xterm256 == 38) {
      return "00afd7";
    }
    if (xterm256 == 39) {
      return "00afff";
    }
    if (xterm256 == 40) {
      return "00d700";
    }
    if (xterm256 == 41) {
      return "00d75f";
    }
    if (xterm256 == 42) {
      return "00d787";
    }
    if (xterm256 == 43) {
      return "00d7af";
    }
    if (xterm256 == 44) {
      return "00d7d7";
    }
    if (xterm256 == 45) {
      return "00d7ff";
    }
    if (xterm256 == 46) {
      return "00ff00";
    }
    if (xterm256 == 47) {
      return "00ff5f";
    }
    if (xterm256 == 48) {
      return "00ff87";
    }
    if (xterm256 == 49) {
      return "00ffaf";
    }
    if (xterm256 == 50) {
      return "00ffd7";
    }
    if (xterm256 == 51) {
      return "00ffff";
    }
    if (xterm256 == 52) {
      return "5f0000";
    }
    if (xterm256 == 53) {
      return "5f005f";
    }
    if (xterm256 == 54) {
      return "5f0087";
    }
    if (xterm256 == 55) {
      return "5f00af";
    }
    if (xterm256 == 56) {
      return "5f00d7";
    }
    if (xterm256 == 57) {
      return "5f00ff";
    }
    if (xterm256 == 58) {
      return "5f5f00";
    }
    if (xterm256 == 59) {
      return "5f5f5f";
    }
    if (xterm256 == 60) {
      return "5f5f87";
    }
    if (xterm256 == 61) {
      return "5f5faf";
    }
    if (xterm256 == 62) {
      return "5f5fd7";
    }
    if (xterm256 == 63) {
      return "5f5fff";
    }
    if (xterm256 == 64) {
      return "5f8700";
    }
    if (xterm256 == 65) {
      return "5f875f";
    }
    if (xterm256 == 66) {
      return "5f8787";
    }
    if (xterm256 == 67) {
      return "5f87af";
    }
    if (xterm256 == 68) {
      return "5f87d7";
    }
    if (xterm256 == 69) {
      return "5f87ff";
    }
    if (xterm256 == 70) {
      return "5faf00";
    }
    if (xterm256 == 71) {
      return "5faf5f";
    }
    if (xterm256 == 72) {
      return "5faf87";
    }
    if (xterm256 == 73) {
      return "5fafaf";
    }
    if (xterm256 == 74) {
      return "5fafd7";
    }
    if (xterm256 == 75) {
      return "5fafff";
    }
    if (xterm256 == 76) {
      return "5fd700";
    }
    if (xterm256 == 77) {
      return "5fd75f";
    }
    if (xterm256 == 78) {
      return "5fd787";
    }
    if (xterm256 == 79) {
      return "5fd7af";
    }
    if (xterm256 == 80) {
      return "5fd7d7";
    }
    if (xterm256 == 81) {
      return "5fd7ff";
    }
    if (xterm256 == 82) {
      return "5fff00";
    }
    if (xterm256 == 83) {
      return "5fff5f";
    }
    if (xterm256 == 84) {
      return "5fff87";
    }
    if (xterm256 == 85) {
      return "5fffaf";
    }
    if (xterm256 == 86) {
      return "5fffd7";
    }
    if (xterm256 == 87) {
      return "5fffff";
    }
    if (xterm256 == 88) {
      return "870000";
    }
    if (xterm256 == 89) {
      return "87005f";
    }
    if (xterm256 == 90) {
      return "870087";
    }
    if (xterm256 == 91) {
      return "8700af";
    }
    if (xterm256 == 92) {
      return "8700d7";
    }
    if (xterm256 == 93) {
      return "8700ff";
    }
    if (xterm256 == 94) {
      return "875f00";
    }
    if (xterm256 == 95) {
      return "875f5f";
    }
    if (xterm256 == 96) {
      return "875f87";
    }
    if (xterm256 == 97) {
      return "875faf";
    }
    if (xterm256 == 98) {
      return "875fd7";
    }
    if (xterm256 == 99) {
      return "875fff";
    }
    if (xterm256 == 100) {
      return "878700";
    }
    if (xterm256 == 101) {
      return "87875f";
    }
    if (xterm256 == 102) {
      return "878787";
    }
    if (xterm256 == 103) {
      return "8787af";
    }
    if (xterm256 == 104) {
      return "8787d7";
    }
    if (xterm256 == 105) {
      return "8787ff";
    }
    if (xterm256 == 106) {
      return "87af00";
    }
    if (xterm256 == 107) {
      return "87af5f";
    }
    if (xterm256 == 108) {
      return "87af87";
    }
    if (xterm256 == 109) {
      return "87afaf";
    }
    if (xterm256 == 110) {
      return "87afd7";
    }
    if (xterm256 == 111) {
      return "87afff";
    }
    if (xterm256 == 112) {
      return "87d700";
    }
    if (xterm256 == 113) {
      return "87d75f";
    }
    if (xterm256 == 114) {
      return "87d787";
    }
    if (xterm256 == 115) {
      return "87d7af";
    }
    if (xterm256 == 116) {
      return "87d7d7";
    }
    if (xterm256 == 117) {
      return "87d7ff";
    }
    if (xterm256 == 118) {
      return "87ff00";
    }
    if (xterm256 == 119) {
      return "87ff5f";
    }
    if (xterm256 == 120) {
      return "87ff87";
    }
    if (xterm256 == 121) {
      return "87ffaf";
    }
    if (xterm256 == 122) {
      return "87ffd7";
    }
    if (xterm256 == 123) {
      return "87ffff";
    }
    if (xterm256 == 124) {
      return "af0000";
    }
    if (xterm256 == 125) {
      return "af005f";
    }
    if (xterm256 == 126) {
      return "af0087";
    }
    if (xterm256 == 127) {
      return "af00af";
    }
    if (xterm256 == 128) {
      return "af00d7";
    }
    if (xterm256 == 129) {
      return "af00ff";
    }
    if (xterm256 == 130) {
      return "af5f00";
    }
    if (xterm256 == 131) {
      return "af5f5f";
    }
    if (xterm256 == 132) {
      return "af5f87";
    }
    if (xterm256 == 133) {
      return "af5faf";
    }
    if (xterm256 == 134) {
      return "af5fd7";
    }
    if (xterm256 == 135) {
      return "af5fff";
    }
    if (xterm256 == 136) {
      return "af8700";
    }
    if (xterm256 == 137) {
      return "af875f";
    }
    if (xterm256 == 138) {
      return "af8787";
    }
    if (xterm256 == 139) {
      return "af87af";
    }
    if (xterm256 == 140) {
      return "af87d7";
    }
    if (xterm256 == 141) {
      return "af87ff";
    }
    if (xterm256 == 142) {
      return "afaf00";
    }
    if (xterm256 == 143) {
      return "afaf5f";
    }
    if (xterm256 == 144) {
      return "afaf87";
    }
    if (xterm256 == 145) {
      return "afafaf";
    }
    if (xterm256 == 146) {
      return "afafd7";
    }
    if (xterm256 == 147) {
      return "afafff";
    }
    if (xterm256 == 148) {
      return "afd700";
    }
    if (xterm256 == 149) {
      return "afd75f";
    }
    if (xterm256 == 150) {
      return "afd787";
    }
    if (xterm256 == 151) {
      return "afd7af";
    }
    if (xterm256 == 152) {
      return "afd7d7";
    }
    if (xterm256 == 153) {
      return "afd7ff";
    }
    if (xterm256 == 154) {
      return "afff00";
    }
    if (xterm256 == 155) {
      return "afff5f";
    }
    if (xterm256 == 156) {
      return "afff87";
    }
    if (xterm256 == 157) {
      return "afffaf";
    }
    if (xterm256 == 158) {
      return "afffd7";
    }
    if (xterm256 == 159) {
      return "afffff";
    }
    if (xterm256 == 160) {
      return "d70000";
    }
    if (xterm256 == 161) {
      return "d7005f";
    }
    if (xterm256 == 162) {
      return "d70087";
    }
    if (xterm256 == 163) {
      return "d700af";
    }
    if (xterm256 == 164) {
      return "d700d7";
    }
    if (xterm256 == 165) {
      return "d700ff";
    }
    if (xterm256 == 166) {
      return "d75f00";
    }
    if (xterm256 == 167) {
      return "d75f5f";
    }
    if (xterm256 == 168) {
      return "d75f87";
    }
    if (xterm256 == 169) {
      return "d75faf";
    }
    if (xterm256 == 170) {
      return "d75fd7";
    }
    if (xterm256 == 171) {
      return "d75fff";
    }
    if (xterm256 == 172) {
      return "d78700";
    }
    if (xterm256 == 173) {
      return "d7875f";
    }
    if (xterm256 == 174) {
      return "d78787";
    }
    if (xterm256 == 175) {
      return "d787af";
    }
    if (xterm256 == 176) {
      return "d787d7";
    }
    if (xterm256 == 177) {
      return "d787ff";
    }
    if (xterm256 == 178) {
      return "d7af00";
    }
    if (xterm256 == 179) {
      return "d7af5f";
    }
    if (xterm256 == 180) {
      return "d7af87";
    }
    if (xterm256 == 181) {
      return "d7afaf";
    }
    if (xterm256 == 182) {
      return "d7afd7";
    }
    if (xterm256 == 183) {
      return "d7afff";
    }
    if (xterm256 == 184) {
      return "d7d700";
    }
    if (xterm256 == 185) {
      return "d7d75f";
    }
    if (xterm256 == 186) {
      return "d7d787";
    }
    if (xterm256 == 187) {
      return "d7d7af";
    }
    if (xterm256 == 188) {
      return "d7d7d7";
    }
    if (xterm256 == 189) {
      return "d7d7ff";
    }
    if (xterm256 == 190) {
      return "d7ff00";
    }
    if (xterm256 == 191) {
      return "d7ff5f";
    }
    if (xterm256 == 192) {
      return "d7ff87";
    }
    if (xterm256 == 193) {
      return "d7ffaf";
    }
    if (xterm256 == 194) {
      return "d7ffd7";
    }
    if (xterm256 == 195) {
      return "d7ffff";
    }
    if (xterm256 == 196) {
      return "ff0000";
    }
    if (xterm256 == 197) {
      return "ff005f";
    }
    if (xterm256 == 198) {
      return "ff0087";
    }
    if (xterm256 == 199) {
      return "ff00af";
    }
    if (xterm256 == 200) {
      return "ff00d7";
    }
    if (xterm256 == 201) {
      return "ff00ff";
    }
    if (xterm256 == 202) {
      return "ff5f00";
    }
    if (xterm256 == 203) {
      return "ff5f5f";
    }
    if (xterm256 == 204) {
      return "ff5f87";
    }
    if (xterm256 == 205) {
      return "ff5faf";
    }
    if (xterm256 == 206) {
      return "ff5fd7";
    }
    if (xterm256 == 207) {
      return "ff5fff";
    }
    if (xterm256 == 208) {
      return "ff8700";
    }
    if (xterm256 == 209) {
      return "ff875f";
    }
    if (xterm256 == 210) {
      return "ff8787";
    }
    if (xterm256 == 211) {
      return "ff87af";
    }
    if (xterm256 == 212) {
      return "ff87d7";
    }
    if (xterm256 == 213) {
      return "ff87ff";
    }
    if (xterm256 == 214) {
      return "ffaf00";
    }
    if (xterm256 == 215) {
      return "ffaf5f";
    }
    if (xterm256 == 216) {
      return "ffaf87";
    }
    if (xterm256 == 217) {
      return "ffafaf";
    }
    if (xterm256 == 218) {
      return "ffafd7";
    }
    if (xterm256 == 219) {
      return "ffafff";
    }
    if (xterm256 == 220) {
      return "ffd700";
    }
    if (xterm256 == 221) {
      return "ffd75f";
    }
    if (xterm256 == 222) {
      return "ffd787";
    }
    if (xterm256 == 223) {
      return "ffd7af";
    }
    if (xterm256 == 224) {
      return "ffd7d7";
    }
    if (xterm256 == 225) {
      return "ffd7ff";
    }
    if (xterm256 == 226) {
      return "ffff00";
    }
    if (xterm256 == 227) {
      return "ffff5f";
    }
    if (xterm256 == 228) {
      return "ffff87";
    }
    if (xterm256 == 229) {
      return "ffffaf";
    }
    if (xterm256 == 230) {
      return "ffffd7";
    }
    if (xterm256 == 231) {
      return "ffffff";
    }

    // Gray-scale range.
    if (xterm256 == 232) {
      return "080808";
    }
    if (xterm256 == 233) {
      return "121212";
    }
    if (xterm256 == 234) {
      return "1c1c1c";
    }
    if (xterm256 == 235) {
      return "262626";
    }
    if (xterm256 == 236) {
      return "303030";
    }
    if (xterm256 == 237) {
      return "3a3a3a";
    }
    if (xterm256 == 238) {
      return "444444";
    }
    if (xterm256 == 239) {
      return "4e4e4e";
    }
    if (xterm256 == 240) {
      return "585858";
    }
    if (xterm256 == 241) {
      return "626262";
    }
    if (xterm256 == 242) {
      return "6c6c6c";
    }
    if (xterm256 == 243) {
      return "767676";
    }
    if (xterm256 == 244) {
      return "808080";
    }
    if (xterm256 == 245) {
      return "8a8a8a";
    }
    if (xterm256 == 246) {
      return "949494";
    }
    if (xterm256 == 247) {
      return "9e9e9e";
    }
    if (xterm256 == 248) {
      return "a8a8a8";
    }
    if (xterm256 == 249) {
      return "b2b2b2";
    }
    if (xterm256 == 250) {
      return "bcbcbc";
    }
    if (xterm256 == 251) {
      return "c6c6c6";
    }
    if (xterm256 == 252) {
      return "d0d0d0";
    }
    if (xterm256 == 253) {
      return "dadada";
    }
    if (xterm256 == 254) {
      return "e4e4e4";
    }
    if (xterm256 == 255) {
      return "eeeeee";
    }

    return null;
  }

  // Map<Integer, String> xterm256toRGB = new HashMap<Integer, String>();

  /**
   * Return the colour name contrasting with the input colour. See Xterm256.colourNameToXterm256()
   * for how input colour can be specified.
   *
   * @throws InvalidColourException
   */
  public static String getContrastColour(String colour) throws InvalidColourException {
    int icol = colourNameToXterm256(colour);
    return contrastColour.get(icol);
  }
}
