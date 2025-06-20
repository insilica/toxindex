import React, { useEffect, useRef, useState } from 'react';

interface TypewriterTextProps {
  text: string;
  typingSpeed?: number;
  restartDelay?: number;
  cursorColor?: string;
  className?: string;
}

const TypewriterText: React.FC<TypewriterTextProps> = ({
  text,
  typingSpeed = 80,
  restartDelay = 5000,
  cursorColor = '#6ee7b7',
  className = ''
}) => {
  const [typedText, setTypedText] = useState("");
  const [isTyping, setIsTyping] = useState(true);
  const typewriterTimeoutRef = useRef<number | null>(null);
  const typewriterIntervalRef = useRef<number | null>(null);

  useEffect(() => {
    function startTypewriter() {
      let i = 0;
      setTypedText("");
      setIsTyping(true);
      if (typewriterIntervalRef.current) clearInterval(typewriterIntervalRef.current);
      if (typewriterTimeoutRef.current) clearTimeout(typewriterTimeoutRef.current);
      
      typewriterIntervalRef.current = window.setInterval(() => {
        setTypedText(text.slice(0, i + 1));
        i++;
        if (i === text.length) {
          if (typewriterIntervalRef.current) clearInterval(typewriterIntervalRef.current);
          setIsTyping(false);
          // Restart after delay
          typewriterTimeoutRef.current = window.setTimeout(() => {
            startTypewriter();
          }, restartDelay);
        }
      }, typingSpeed);
    }

    startTypewriter();
    return () => {
      if (typewriterIntervalRef.current) clearInterval(typewriterIntervalRef.current);
      if (typewriterTimeoutRef.current) clearTimeout(typewriterTimeoutRef.current);
    };
  }, [text, typingSpeed, restartDelay]);

  return (
    <span className={className}>
      {typedText}
      {isTyping && (
        <span 
          style={{ 
            color: cursorColor, 
            fontWeight: 'bold', 
            marginLeft: 2, 
            animation: 'blink 1s steps(1) infinite' 
          }}
        >
          |
        </span>
      )}
      <style>{`
        @keyframes blink {
          0%, 100% { opacity: 1; }
          50% { opacity: 0; }
        }
      `}</style>
    </span>
  );
};

export default TypewriterText; 