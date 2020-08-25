package utils;

import static org.junit.Assert.*;

import java.util.List;

import org.junit.Test;

public class TokenizerTest {

	@Test
	public void canTokenizeEmptyString() {
		
		List<String> tk= new Tokenizer("").tokenize();
		assertEquals(0, tk.size());

		tk= new Tokenizer("  ").tokenize();
		assertEquals(0, tk.size());
	}
	
	@Test
	public void canTokenizeUnquotedStrings() {

		String cmd= "do   -c AA stuff";
		List<String> tk= new Tokenizer(cmd).tokenize();;
		assertEquals(4, tk.size());
		assertEquals("AA", tk.get(2));
		
		cmd= "A   -c foo Z";
		tk= new Tokenizer(cmd).tokenize();;
		assertEquals(4, tk.size());
		assertEquals("A", tk.get(0));
		assertEquals("Z", tk.get(3));
		
		cmd= "  do  ";
		tk= new Tokenizer(cmd).tokenize();;
		assertEquals(1, tk.size());
		assertEquals("do", tk.get(0));
		
		cmd= "  do   -c this    stuff  ";
		tk= new Tokenizer(cmd).tokenize();
		assertEquals(4, tk.size());
		assertEquals("do", tk.get(0));
		assertEquals("stuff", tk.get(3));
	}
	
	@Test
	public void canTokenizeQuotes() {
		String cmd= "A 'B' C";
		List<String> tk= new Tokenizer(cmd).tokenize();;
		assertEquals(3, tk.size());
		assertEquals("B", tk.get(1));
		
		cmd= "do   -c 'this    stuff  '";
		tk= new Tokenizer(cmd).tokenize();;
		assertEquals("do", tk.get(0));
		assertEquals("-c", tk.get(1));
		assertEquals("this    stuff  ", tk.get(2));

		cmd= "do   -c 'this    stuff  ' ' foo bar'";
		tk= new Tokenizer(cmd).tokenize();
		assertEquals(4, tk.size());
		assertEquals("do", tk.get(0));
		assertEquals("this    stuff  ", tk.get(2));
		assertEquals(" foo bar", tk.get(3));
		
		cmd= "do   -c \"this    stuff  \"";
		tk= new Tokenizer(cmd).tokenize();;
		assertEquals("do", tk.get(0));
		assertEquals("-c", tk.get(1));
		assertEquals("this    stuff  ", tk.get(2));
		
		cmd= "do   -c '''this    stuff  '''";
		tk= new Tokenizer(cmd).tokenize();;
		assertEquals("do", tk.get(0));
		assertEquals("-c", tk.get(1));
		assertEquals("this    stuff  ", tk.get(2));
		
		cmd= "do   -c \"\"\"this    stuff  \"\"\"";
		tk= new Tokenizer(cmd).tokenize();;
		assertEquals("do", tk.get(0));
		assertEquals("-c", tk.get(1));
		assertEquals("this    stuff  ", tk.get(2));
	}
	
	@Test
	public void canTokenizeStringWithEmptyToken() {
		
		String cmd= "do   -c ' ' stuff ''' '''";
		List<String> tk= new Tokenizer(cmd).tokenize();;
		assertEquals(5, tk.size());
		assertEquals(" ", tk.get(2));
		assertEquals(" ", tk.get(4));
		
		cmd= "do   -c '' stuff ''";
		tk= new Tokenizer(cmd).tokenize();
		assertEquals(5, tk.size());
		assertEquals("", tk.get(2));
		assertEquals("", tk.get(4));
		
		cmd= "do   -c '''''' stuff \"\"\"\"\"\"";
		tk= new Tokenizer(cmd).tokenize();
		assertEquals(5, tk.size());
		assertEquals("", tk.get(2));
		assertEquals("", tk.get(4));
	}

	@Test
	public void canTokenizeNestedQuotes() {
		String cmd= "do 'foo \" \tbar'";
		List<String> tk= new Tokenizer(cmd).tokenize();
		assertEquals(2, tk.size());
		assertEquals("foo \" \tbar", tk.get(1));
		
		cmd= "do '''foo \\' bar'''";
		tk= new Tokenizer(cmd).tokenize();
		System.err.println(tk);
	}
	@Test
	public void separatorBetweenQuotedTokensIsOptional() {
		// Is this behaviour desirable?
		String cmd= "'foo' 'bar' 'baz'";
		String cmd2= "'foo''bar''baz'";
		List<String> tk= new Tokenizer(cmd).tokenize();
		List<String> tk2= new Tokenizer(cmd2).tokenize();
		assertEquals(tk, tk2);
		
	}
	
	@Test
	public void canEscape() {
		String cmd= "do 'foo\\\\\\'bar'";
		List<String> tk= new Tokenizer(cmd).tokenize();;
		assertEquals(2, tk.size());
		assertEquals("foo\\\\\\'bar", tk.get(1));
		
		cmd= "do 'foo\\'bar'";
		tk= new Tokenizer(cmd).tokenize();;
		assertEquals(2, tk.size());
		assertEquals("foo\\'bar", tk.get(1));
		
		cmd= "do 'foo\\\\' bar";
		tk= new Tokenizer(cmd).tokenize();;
		assertEquals(3, tk.size());
		assertEquals("foo\\\\", tk.get(1));
		assertEquals("bar", tk.get(2));
		
	}
	
}
