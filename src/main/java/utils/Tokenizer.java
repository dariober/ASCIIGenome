package utils;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.apache.commons.lang3.StringUtils;

public class Tokenizer {

	private char SEP= ' '; 
	private Map<String, Boolean> quoteChars=  new HashMap<String, Boolean>();
	private int cursor= 0;
	private String string;
	private static final String S3Q= "'''";
	private static final String D3Q= "\"\"\"";
	private static final String SQ= "'";
	private static final String DQ= "\"";
	
	public Tokenizer(String string){
		this.string= string.trim();
		quoteChars.put(SQ, false);
		quoteChars.put(DQ, false);
		quoteChars.put(S3Q, false);
		quoteChars.put(D3Q, false);
	}
	
	public List<String> tokenize() {

		List<String> tokens= new ArrayList<String>();

		StringBuilder token= null;
		while(this.cursor < this.string.length()) {
			Character c= this.string.charAt(this.cursor);
			if(token == null) {
				// The char we just read is at the beginning of a substring (token)
				// Is this char a quote? If so, record the fact that we are starting a quoted string.
				token= new StringBuilder();
				boolean isQuoted= this.openQuote();
				if(isQuoted) {
					continue;
				}
			} 
			else {
				// We are inside a growing token. Check the char we just read matches an open quote, 
				// if any. If so, terminate the token and move cursor.
				boolean isQuoted= this.closeQuote();
				if(isQuoted) {
					tokens.add(token.toString());
					token= null;
					continue;
				}
			}
			if(this.isSeparator()) {
				tokens.add(token.toString());
				token= null;
			} 
			else {
				token.append(c);
				this.cursor += 1;
			}
		}
		if(token != null && token.length() > 0) {
			tokens.add(token.toString());
		}
		return tokens;
	}

	private boolean closeQuote() {
		String c= String.valueOf(this.string.charAt(this.cursor));
		String quote= "";
		for(String x : this.quoteChars.keySet()) {
			if(this.quoteChars.get(x)) {
				quote= x;
				break;
			}
		}
		if(quote.isEmpty()) {
			// We are not inside an open quote
			return false;
		}
		
		if( ! c.equals(DQ) && ! c.equals(SQ) ) {
			// The char we are inspecting is not a quote
			return false;
		}
		
		// Check whether the quoting character is escaped 
		String pre= this.string.substring(0, this.cursor);
		String pre2= StringUtils.stripEnd(pre, "\\");
		if((pre.length() - pre2.length()) % 2 == 1) {
			// There is an odd number of escapes preceding this quote char (e.g. 1 backslash)
			// Therefore the quoting is not a closing mark
			return false;
		}
		
		String sub= this.string.substring(this.cursor);
		if(sub.startsWith(quote)) {
			this.quoteChars.put(quote, false);
			// We advance the cursor with the length of the quote itself and possible separators after it.
			// so we place the cursor at the start of the next token
			String strip= sub.replaceAll("^" + quote + this.SEP + "*", "");
			this.cursor += sub.length() - strip.length();
			return true;			
		} else {
			return false;
		}
	}
	
	private boolean  openQuote() {
		String c= String.valueOf(this.string.charAt(this.cursor));
		if( ! c.equals(DQ) && ! c.equals(SQ) ) {
			// The substring is not quoted
			return false;
		}
		// What sort of quoting we have single? double? Triple?
		// Move the cursor accordingly
		String quote= "";
		String sub= this.string.substring(this.cursor);
		if(sub.startsWith(S3Q)) {
			quote= S3Q;
		} else if(sub.startsWith(D3Q)) {
			quote= D3Q;
		} else if(sub.startsWith(DQ)) {
			quote= DQ;
		} else if(sub.startsWith(SQ)) {
			quote= SQ;
		}
		this.quoteChars.put(quote, true);
		// A sanity check that one and only one quote is open
		int n= 0;
		for(boolean x : this.quoteChars.values()) {
			if(x) {
				n += 1;
			}
		}
		if(n != 1) {
			throw new RuntimeException();
		}
		this.cursor += quote.length();
		return true;
	}
	
	private boolean isSeparator() {
		for(boolean x : this.quoteChars.values()) {
			if(x) {
				// We are inside a quoted string so ignore spaces.
				return false;
			}
		}
		if(this.string.charAt(this.cursor) != this.SEP) {
			return false;
		} 
		// Current position IS a delimiter
		for(int i= this.cursor; i < this.string.length(); i++) {
			// Advance cursor until the next non-space char
			if(this.string.charAt(this.cursor) != this.SEP) {
				break;
			}
			this.cursor= i;
		}
		return true;
	}
	
}
