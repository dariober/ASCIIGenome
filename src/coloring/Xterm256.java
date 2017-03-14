package coloring;

import java.awt.Color;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.Map;

import org.apache.commons.lang3.StringUtils;

import exceptions.InvalidColourException;

public class Xterm256 {

	public static final Color xterm256ToColor(int xterm256) throws InvalidColourException{
		
		String rgb= xterm256toRGB(xterm256);
		
		int r= Integer.parseInt(rgb.substring(0, 2), 16);
		int g= Integer.parseInt(rgb.substring(2, 4), 16);
		int b= Integer.parseInt(rgb.substring(4, 6), 16);
		
		float[] hsb= Color.RGBtoHSB(r, g, b, null);
		return Color.getHSBColor(hsb[0], hsb[1], hsb[2]);
		
	}
	
	public static final LinkedHashMap<String, Integer> mapColorNameToXterm256(){
		
		final LinkedHashMap<String, Integer> map= new LinkedHashMap<String, Integer>();
		
		// See http://jonasjacek.github.io/colors/
		map.put("black", 0);
		map.put("maroon", 1);
		map.put("green", 2);
		map.put("olive", 3);
		map.put("navy", 4);
		map.put("purple", 5);
		map.put("teal", 6);
		map.put("silver", 7);
		map.put("grey", 8);
		map.put("red", 9);
		map.put("lime", 10);
		map.put("yellow", 11);
		map.put("blue", 12);
		map.put("fuchsia", 13);
		map.put("aqua", 14);
		map.put("white", 231); // NB: white is "15" but 15 shows up as grey-ish. So use 231 (grey100) instead 
		map.put("grey0", 16);
		map.put("navyblue", 17);
		map.put("darkblue", 18);
		map.put("blue3", 19);
		map.put("blue3", 20);
		map.put("blue1", 21);
		map.put("darkgreen", 22);
		map.put("deepskyblue4", 23);
		map.put("deepskyblue4", 24);
		map.put("deepskyblue4", 25);
		map.put("dodgerblue3", 26);
		map.put("dodgerblue2", 27);
		map.put("green4", 28);
		map.put("springgreen4", 29);
		map.put("turquoise4", 30);
		map.put("deepskyblue3", 31);
		map.put("deepskyblue3", 32);
		map.put("dodgerblue1", 33);
		map.put("green3", 34);
		map.put("springgreen3", 35);
		map.put("darkcyan", 36);
		map.put("lightseagreen", 37);
		map.put("deepskyblue2", 38);
		map.put("deepskyblue1", 39);
		map.put("green3", 40);
		map.put("springgreen3", 41);
		map.put("springgreen2", 42);
		map.put("cyan3", 43);
		map.put("darkturquoise", 44);
		map.put("turquoise2", 45);
		map.put("green1", 46);
		map.put("springgreen2", 47);
		map.put("springgreen1", 48);
		map.put("mediumspringgreen", 49);
		map.put("cyan2", 50);
		map.put("cyan1", 51);
		map.put("darkred", 52);
		map.put("deeppink4", 53);
		map.put("purple4", 54);
		map.put("purple4", 55);
		map.put("purple3", 56);
		map.put("blueviolet", 57);
		map.put("orange4", 58);
		map.put("grey37", 59);
		map.put("mediumpurple4", 60);
		map.put("slateblue3", 61);
		map.put("slateblue3", 62);
		map.put("royalblue1", 63);
		map.put("chartreuse4", 64);
		map.put("darkseagreen4", 65);
		map.put("paleturquoise4", 66);
		map.put("steelblue", 67);
		map.put("steelblue3", 68);
		map.put("cornflowerblue", 69);
		map.put("chartreuse3", 70);
		map.put("darkseagreen4", 71);
		map.put("cadetblue", 72);
		map.put("cadetblue", 73);
		map.put("skyblue3", 74);
		map.put("steelblue1", 75);
		map.put("chartreuse3", 76);
		map.put("palegreen3", 77);
		map.put("seagreen3", 78);
		map.put("aquamarine3", 79);
		map.put("mediumturquoise", 80);
		map.put("steelblue1", 81);
		map.put("chartreuse2", 82);
		map.put("seagreen2", 83);
		map.put("seagreen1", 84);
		map.put("seagreen1", 85);
		map.put("aquamarine1", 86);
		map.put("darkslategray2", 87);
		map.put("darkred", 88);
		map.put("deeppink4", 89);
		map.put("darkmagenta", 90);
		map.put("darkmagenta", 91);
		map.put("darkviolet", 92);
		map.put("purple", 93);
		map.put("orange4", 94);
		map.put("lightpink4", 95);
		map.put("plum4", 96);
		map.put("mediumpurple3", 97);
		map.put("mediumpurple3", 98);
		map.put("slateblue1", 99);
		map.put("yellow4", 100);
		map.put("wheat4", 101);
		map.put("grey53", 102);
		map.put("lightslategrey", 103);
		map.put("mediumpurple", 104);
		map.put("lightslateblue", 105);
		map.put("yellow4", 106);
		map.put("darkolivegreen3", 107);
		map.put("darkseagreen", 108);
		map.put("lightskyblue3", 109);
		map.put("lightskyblue3", 110);
		map.put("skyblue2", 111);
		map.put("chartreuse2", 112);
		map.put("darkolivegreen3", 113);
		map.put("palegreen3", 114);
		map.put("darkseagreen3", 115);
		map.put("darkslategray3", 116);
		map.put("skyblue1", 117);
		map.put("chartreuse1", 118);
		map.put("lightgreen", 119);
		map.put("lightgreen", 120);
		map.put("palegreen1", 121);
		map.put("aquamarine1", 122);
		map.put("darkslategray1", 123);
		map.put("red3", 124);
		map.put("deeppink4", 125);
		map.put("mediumvioletred", 126);
		map.put("magenta3", 127);
		map.put("darkviolet", 128);
		map.put("purple", 129);
		map.put("darkorange3", 130);
		map.put("indianred", 131);
		map.put("hotpink3", 132);
		map.put("mediumorchid3", 133);
		map.put("mediumorchid", 134);
		map.put("mediumpurple2", 135);
		map.put("darkgoldenrod", 136);
		map.put("lightsalmon3", 137);
		map.put("rosybrown", 138);
		map.put("grey63", 139);
		map.put("mediumpurple2", 140);
		map.put("mediumpurple1", 141);
		map.put("gold3", 142);
		map.put("darkkhaki", 143);
		map.put("navajowhite3", 144);
		map.put("grey69", 145);
		map.put("lightsteelblue3", 146);
		map.put("lightsteelblue", 147);
		map.put("yellow3", 148);
		map.put("darkolivegreen3", 149);
		map.put("darkseagreen3", 150);
		map.put("darkseagreen2", 151);
		map.put("lightcyan3", 152);
		map.put("lightskyblue1", 153);
		map.put("greenyellow", 154);
		map.put("darkolivegreen2", 155);
		map.put("palegreen1", 156);
		map.put("darkseagreen2", 157);
		map.put("darkseagreen1", 158);
		map.put("paleturquoise1", 159);
		map.put("red3", 160);
		map.put("deeppink3", 161);
		map.put("deeppink3", 162);
		map.put("magenta3", 163);
		map.put("magenta3", 164);
		map.put("magenta2", 165);
		map.put("darkorange3", 166);
		map.put("indianred", 167);
		map.put("hotpink3", 168);
		map.put("hotpink2", 169);
		map.put("orchid", 170);
		map.put("mediumorchid1", 171);
		map.put("orange3", 172);
		map.put("lightsalmon3", 173);
		map.put("lightpink3", 174);
		map.put("pink3", 175);
		map.put("plum3", 176);
		map.put("violet", 177);
		map.put("gold3", 178);
		map.put("lightgoldenrod3", 179);
		map.put("tan", 180);
		map.put("mistyrose3", 181);
		map.put("thistle3", 182);
		map.put("plum2", 183);
		map.put("yellow3", 184);
		map.put("khaki3", 185);
		map.put("lightgoldenrod2", 186);
		map.put("lightyellow3", 187);
		map.put("grey84", 188);
		map.put("lightsteelblue1", 189);
		map.put("yellow2", 190);
		map.put("darkolivegreen1", 191);
		map.put("darkolivegreen1", 192);
		map.put("darkseagreen1", 193);
		map.put("honeydew2", 194);
		map.put("lightcyan1", 195);
		map.put("red1", 196);
		map.put("deeppink2", 197);
		map.put("deeppink1", 198);
		map.put("deeppink1", 199);
		map.put("magenta2", 200);
		map.put("magenta1", 201);
		map.put("orangered1", 202);
		map.put("indianred1", 203);
		map.put("indianred1", 204);
		map.put("hotpink", 205);
		map.put("hotpink", 206);
		map.put("mediumorchid1", 207);
		map.put("darkorange", 208);
		map.put("salmon1", 209);
		map.put("lightcoral", 210);
		map.put("palevioletred1", 211);
		map.put("orchid2", 212);
		map.put("orchid1", 213);
		map.put("orange1", 214);
		map.put("sandybrown", 215);
		map.put("lightsalmon1", 216);
		map.put("lightpink1", 217);
		map.put("pink1", 218);
		map.put("plum1", 219);
		map.put("gold1", 220);
		map.put("lightgoldenrod2", 221);
		map.put("lightgoldenrod2", 222);
		map.put("navajowhite1", 223);
		map.put("mistyrose1", 224);
		map.put("thistle1", 225);
		map.put("yellow1", 226);
		map.put("lightgoldenrod1", 227);
		map.put("khaki1", 228);
		map.put("wheat1", 229);
		map.put("cornsilk1", 230);
		map.put("grey100", 231);
		map.put("grey3", 232);
		map.put("grey7", 233);
		map.put("grey11", 234);
		map.put("grey15", 235);
		map.put("grey19", 236);
		map.put("grey23", 237);
		map.put("grey27", 238);
		map.put("grey30", 239);
		map.put("grey35", 240);
		map.put("grey39", 241);
		map.put("grey42", 242);
		map.put("grey46", 243);
		map.put("grey50", 244);
		map.put("grey54", 245);
		map.put("grey58", 246);
		map.put("grey62", 247);
		map.put("grey66", 248);
		map.put("grey70", 249);
		map.put("grey74", 250);
		map.put("grey78", 251);
		map.put("grey82", 252);
		map.put("grey85", 253);
		map.put("grey89", 254);
		map.put("grey93", 255); 

		return map;
	}
	
	public static int colorNameToXterm256(String colorName) throws InvalidColourException {
			
		try{ // See if colour name is already an int between 0 and 255
			int xterm= Integer.parseInt(colorName);
			if(xterm >= 0 && xterm  <= 255){
				return xterm;
			}
		} catch (NumberFormatException e){
			//
		}
		
		colorName= colorName.toLowerCase();
		
        if(mapColorNameToXterm256().containsKey(colorName)) {
        	return mapColorNameToXterm256().get(colorName);
        } else {
        	// Try approximate matching
        	for(String x : mapColorNameToXterm256().keySet()){
        		if(x.startsWith(colorName)){
        			return mapColorNameToXterm256().get(x);
        		}
        	}
			System.err.println("Unrecognized color: " + colorName);
			throw new InvalidColourException();
		}		
	}
	
	public static String colorShowForTerminal() throws InvalidColourException{
		
		int maxLen= 0;
		for(String x : mapColorNameToXterm256().keySet()){
			if(x.length() > maxLen){
				maxLen= x.length(); 
			}
		}
		
		StringBuilder sb= new StringBuilder();
		int i= 0;
		for(String x : mapColorNameToXterm256().keySet()){
			i++;
			int xterm= mapColorNameToXterm256().get(x);
			int spacer= maxLen - x.length();
			sb.append(xterm + ": \033[38;5;" + xterm + "m" + x + 
					"\033[38;5;" + Config.getColor(ConfigKey.foreground) + 
					";48;5;" + Config.getColor(ConfigKey.background) + "m");
			
			if(i == 3){ // Arrange colors in this many columns
				sb.append("\n");
				i= 0;
			} else {
				sb.append(StringUtils.repeat(" ", spacer + 1));
			}
		}
		return sb.toString();
	}
	
	private static final String xterm256toRGB(int xterm256) throws InvalidColourException{
		
		if(xterm256 < 0 || xterm256 > 255){
			throw new InvalidColourException();
		}
		
		// From https://gist.github.com/MicahElliott/719710
		
		// Primary 3-bit (8 colors). Unique representation!
		if( xterm256 == 0) { return "000000";}
		if( xterm256 == 1) { return "800000";}
		if( xterm256 == 2) { return "008000";}
		if( xterm256 == 3) { return "808000";}
		if( xterm256 == 4) { return "000080";}
		if( xterm256 == 5) { return "800080";}
		if( xterm256 == 6) { return "008080";}
		if( xterm256 == 7) { return "c0c0c0";}

		// Equivalent "bright" versions of original 8 colors.
		if( xterm256 == 8) { return "808080";}
		if( xterm256 == 9) { return "ff0000";}
		if( xterm256 == 10) { return "00ff00";}
		if( xterm256 == 11) { return "ffff00";}
		if( xterm256 == 12) { return "0000ff";}
		if( xterm256 == 13) { return "ff00ff";}
		if( xterm256 == 14) { return "00ffff";}
		if( xterm256 == 15) { return "ffffff";}

		// Strictly ascending.
		if( xterm256 == 16) { return "000000";}
		if( xterm256 == 17) { return "00005f";}
		if( xterm256 == 18) { return "000087";}
		if( xterm256 == 19) { return "0000af";}
		if( xterm256 == 20) { return "0000d7";}
		if( xterm256 == 21) { return "0000ff";}
		if( xterm256 == 22) { return "005f00";}
		if( xterm256 == 23) { return "005f5f";}
		if( xterm256 == 24) { return "005f87";}
		if( xterm256 == 25) { return "005faf";}
		if( xterm256 == 26) { return "005fd7";}
		if( xterm256 == 27) { return "005fff";}
		if( xterm256 == 28) { return "008700";}
		if( xterm256 == 29) { return "00875f";}
		if( xterm256 == 30) { return "008787";}
		if( xterm256 == 31) { return "0087af";}
		if( xterm256 == 32) { return "0087d7";}
		if( xterm256 == 33) { return "0087ff";}
		if( xterm256 == 34) { return "00af00";}
		if( xterm256 == 35) { return "00af5f";}
		if( xterm256 == 36) { return "00af87";}
		if( xterm256 == 37) { return "00afaf";}
		if( xterm256 == 38) { return "00afd7";}
		if( xterm256 == 39) { return "00afff";}
		if( xterm256 == 40) { return "00d700";}
		if( xterm256 == 41) { return "00d75f";}
		if( xterm256 == 42) { return "00d787";}
		if( xterm256 == 43) { return "00d7af";}
		if( xterm256 == 44) { return "00d7d7";}
		if( xterm256 == 45) { return "00d7ff";}
		if( xterm256 == 46) { return "00ff00";}
		if( xterm256 == 47) { return "00ff5f";}
		if( xterm256 == 48) { return "00ff87";}
		if( xterm256 == 49) { return "00ffaf";}
		if( xterm256 == 50) { return "00ffd7";}
		if( xterm256 == 51) { return "00ffff";}
		if( xterm256 == 52) { return "5f0000";}
		if( xterm256 == 53) { return "5f005f";}
		if( xterm256 == 54) { return "5f0087";}
		if( xterm256 == 55) { return "5f00af";}
		if( xterm256 == 56) { return "5f00d7";}
		if( xterm256 == 57) { return "5f00ff";}
		if( xterm256 == 58) { return "5f5f00";}
		if( xterm256 == 59) { return "5f5f5f";}
		if( xterm256 == 60) { return "5f5f87";}
		if( xterm256 == 61) { return "5f5faf";}
		if( xterm256 == 62) { return "5f5fd7";}
		if( xterm256 == 63) { return "5f5fff";}
		if( xterm256 == 64) { return "5f8700";}
		if( xterm256 == 65) { return "5f875f";}
		if( xterm256 == 66) { return "5f8787";}
		if( xterm256 == 67) { return "5f87af";}
		if( xterm256 == 68) { return "5f87d7";}
		if( xterm256 == 69) { return "5f87ff";}
		if( xterm256 == 70) { return "5faf00";}
		if( xterm256 == 71) { return "5faf5f";}
		if( xterm256 == 72) { return "5faf87";}
		if( xterm256 == 73) { return "5fafaf";}
		if( xterm256 == 74) { return "5fafd7";}
		if( xterm256 == 75) { return "5fafff";}
		if( xterm256 == 76) { return "5fd700";}
		if( xterm256 == 77) { return "5fd75f";}
		if( xterm256 == 78) { return "5fd787";}
		if( xterm256 == 79) { return "5fd7af";}
		if( xterm256 == 80) { return "5fd7d7";}
		if( xterm256 == 81) { return "5fd7ff";}
		if( xterm256 == 82) { return "5fff00";}
		if( xterm256 == 83) { return "5fff5f";}
		if( xterm256 == 84) { return "5fff87";}
		if( xterm256 == 85) { return "5fffaf";}
		if( xterm256 == 86) { return "5fffd7";}
		if( xterm256 == 87) { return "5fffff";}
		if( xterm256 == 88) { return "870000";}
		if( xterm256 == 89) { return "87005f";}
		if( xterm256 == 90) { return "870087";}
		if( xterm256 == 91) { return "8700af";}
		if( xterm256 == 92) { return "8700d7";}
		if( xterm256 == 93) { return "8700ff";}
		if( xterm256 == 94) { return "875f00";}
		if( xterm256 == 95) { return "875f5f";}
		if( xterm256 == 96) { return "875f87";}
		if( xterm256 == 97) { return "875faf";}
		if( xterm256 == 98) { return "875fd7";}
		if( xterm256 == 99) { return "875fff";}
		if( xterm256 == 100) { return "878700";}
		if( xterm256 == 101) { return "87875f";}
		if( xterm256 == 102) { return "878787";}
		if( xterm256 == 103) { return "8787af";}
		if( xterm256 == 104) { return "8787d7";}
		if( xterm256 == 105) { return "8787ff";}
		if( xterm256 == 106) { return "87af00";}
		if( xterm256 == 107) { return "87af5f";}
		if( xterm256 == 108) { return "87af87";}
		if( xterm256 == 109) { return "87afaf";}
		if( xterm256 == 110) { return "87afd7";}
		if( xterm256 == 111) { return "87afff";}
		if( xterm256 == 112) { return "87d700";}
		if( xterm256 == 113) { return "87d75f";}
		if( xterm256 == 114) { return "87d787";}
		if( xterm256 == 115) { return "87d7af";}
		if( xterm256 == 116) { return "87d7d7";}
		if( xterm256 == 117) { return "87d7ff";}
		if( xterm256 == 118) { return "87ff00";}
		if( xterm256 == 119) { return "87ff5f";}
		if( xterm256 == 120) { return "87ff87";}
		if( xterm256 == 121) { return "87ffaf";}
		if( xterm256 == 122) { return "87ffd7";}
		if( xterm256 == 123) { return "87ffff";}
		if( xterm256 == 124) { return "af0000";}
		if( xterm256 == 125) { return "af005f";}
		if( xterm256 == 126) { return "af0087";}
		if( xterm256 == 127) { return "af00af";}
		if( xterm256 == 128) { return "af00d7";}
		if( xterm256 == 129) { return "af00ff";}
		if( xterm256 == 130) { return "af5f00";}
		if( xterm256 == 131) { return "af5f5f";}
		if( xterm256 == 132) { return "af5f87";}
		if( xterm256 == 133) { return "af5faf";}
		if( xterm256 == 134) { return "af5fd7";}
		if( xterm256 == 135) { return "af5fff";}
		if( xterm256 == 136) { return "af8700";}
		if( xterm256 == 137) { return "af875f";}
		if( xterm256 == 138) { return "af8787";}
		if( xterm256 == 139) { return "af87af";}
		if( xterm256 == 140) { return "af87d7";}
		if( xterm256 == 141) { return "af87ff";}
		if( xterm256 == 142) { return "afaf00";}
		if( xterm256 == 143) { return "afaf5f";}
		if( xterm256 == 144) { return "afaf87";}
		if( xterm256 == 145) { return "afafaf";}
		if( xterm256 == 146) { return "afafd7";}
		if( xterm256 == 147) { return "afafff";}
		if( xterm256 == 148) { return "afd700";}
		if( xterm256 == 149) { return "afd75f";}
		if( xterm256 == 150) { return "afd787";}
		if( xterm256 == 151) { return "afd7af";}
		if( xterm256 == 152) { return "afd7d7";}
		if( xterm256 == 153) { return "afd7ff";}
		if( xterm256 == 154) { return "afff00";}
		if( xterm256 == 155) { return "afff5f";}
		if( xterm256 == 156) { return "afff87";}
		if( xterm256 == 157) { return "afffaf";}
		if( xterm256 == 158) { return "afffd7";}
		if( xterm256 == 159) { return "afffff";}
		if( xterm256 == 160) { return "d70000";}
		if( xterm256 == 161) { return "d7005f";}
		if( xterm256 == 162) { return "d70087";}
		if( xterm256 == 163) { return "d700af";}
		if( xterm256 == 164) { return "d700d7";}
		if( xterm256 == 165) { return "d700ff";}
		if( xterm256 == 166) { return "d75f00";}
		if( xterm256 == 167) { return "d75f5f";}
		if( xterm256 == 168) { return "d75f87";}
		if( xterm256 == 169) { return "d75faf";}
		if( xterm256 == 170) { return "d75fd7";}
		if( xterm256 == 171) { return "d75fff";}
		if( xterm256 == 172) { return "d78700";}
		if( xterm256 == 173) { return "d7875f";}
		if( xterm256 == 174) { return "d78787";}
		if( xterm256 == 175) { return "d787af";}
		if( xterm256 == 176) { return "d787d7";}
		if( xterm256 == 177) { return "d787ff";}
		if( xterm256 == 178) { return "d7af00";}
		if( xterm256 == 179) { return "d7af5f";}
		if( xterm256 == 180) { return "d7af87";}
		if( xterm256 == 181) { return "d7afaf";}
		if( xterm256 == 182) { return "d7afd7";}
		if( xterm256 == 183) { return "d7afff";}
		if( xterm256 == 184) { return "d7d700";}
		if( xterm256 == 185) { return "d7d75f";}
		if( xterm256 == 186) { return "d7d787";}
		if( xterm256 == 187) { return "d7d7af";}
		if( xterm256 == 188) { return "d7d7d7";}
		if( xterm256 == 189) { return "d7d7ff";}
		if( xterm256 == 190) { return "d7ff00";}
		if( xterm256 == 191) { return "d7ff5f";}
		if( xterm256 == 192) { return "d7ff87";}
		if( xterm256 == 193) { return "d7ffaf";}
		if( xterm256 == 194) { return "d7ffd7";}
		if( xterm256 == 195) { return "d7ffff";}
		if( xterm256 == 196) { return "ff0000";}
		if( xterm256 == 197) { return "ff005f";}
		if( xterm256 == 198) { return "ff0087";}
		if( xterm256 == 199) { return "ff00af";}
		if( xterm256 == 200) { return "ff00d7";}
		if( xterm256 == 201) { return "ff00ff";}
		if( xterm256 == 202) { return "ff5f00";}
		if( xterm256 == 203) { return "ff5f5f";}
		if( xterm256 == 204) { return "ff5f87";}
		if( xterm256 == 205) { return "ff5faf";}
		if( xterm256 == 206) { return "ff5fd7";}
		if( xterm256 == 207) { return "ff5fff";}
		if( xterm256 == 208) { return "ff8700";}
		if( xterm256 == 209) { return "ff875f";}
		if( xterm256 == 210) { return "ff8787";}
		if( xterm256 == 211) { return "ff87af";}
		if( xterm256 == 212) { return "ff87d7";}
		if( xterm256 == 213) { return "ff87ff";}
		if( xterm256 == 214) { return "ffaf00";}
		if( xterm256 == 215) { return "ffaf5f";}
		if( xterm256 == 216) { return "ffaf87";}
		if( xterm256 == 217) { return "ffafaf";}
		if( xterm256 == 218) { return "ffafd7";}
		if( xterm256 == 219) { return "ffafff";}
		if( xterm256 == 220) { return "ffd700";}
		if( xterm256 == 221) { return "ffd75f";}
		if( xterm256 == 222) { return "ffd787";}
		if( xterm256 == 223) { return "ffd7af";}
		if( xterm256 == 224) { return "ffd7d7";}
		if( xterm256 == 225) { return "ffd7ff";}
		if( xterm256 == 226) { return "ffff00";}
		if( xterm256 == 227) { return "ffff5f";}
		if( xterm256 == 228) { return "ffff87";}
		if( xterm256 == 229) { return "ffffaf";}
		if( xterm256 == 230) { return "ffffd7";}
		if( xterm256 == 231) { return "ffffff";}

		// Gray-scale range.
		if( xterm256 == 232) { return "080808";}
		if( xterm256 == 233) { return "121212";}
		if( xterm256 == 234) { return "1c1c1c";}
		if( xterm256 == 235) { return "262626";}
		if( xterm256 == 236) { return "303030";}
		if( xterm256 == 237) { return "3a3a3a";}
		if( xterm256 == 238) { return "444444";}
		if( xterm256 == 239) { return "4e4e4e";}
		if( xterm256 == 240) { return "585858";}
		if( xterm256 == 241) { return "626262";}
		if( xterm256 == 242) { return "6c6c6c";}
		if( xterm256 == 243) { return "767676";}
		if( xterm256 == 244) { return "808080";}
		if( xterm256 == 245) { return "8a8a8a";}
		if( xterm256 == 246) { return "949494";}
		if( xterm256 == 247) { return "9e9e9e";}
		if( xterm256 == 248) { return "a8a8a8";}
		if( xterm256 == 249) { return "b2b2b2";}
		if( xterm256 == 250) { return "bcbcbc";}
		if( xterm256 == 251) { return "c6c6c6";}
		if( xterm256 == 252) { return "d0d0d0";}
		if( xterm256 == 253) { return "dadada";}
		if( xterm256 == 254) { return "e4e4e4";}
		if( xterm256 == 255) { return "eeeeee";}
		
		return null;
		
	}
	
	Map<Integer, String> xterm256toRGB = new HashMap<Integer, String>();

}
