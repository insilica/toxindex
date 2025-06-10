import React, { useEffect, useState } from "react";
import { useLocation } from "react-router-dom";
import ReactMarkdown from "react-markdown";

const TermsPrivacy: React.FC = () => {
  const location = useLocation();
  const isTerms = location.pathname.includes("terms-of-use");
  const [content, setContent] = useState("");

  useEffect(() => {
    const file = isTerms ? "/policies/terms-of-use.md" : "/policies/privacy-policy.md";
    fetch(file)
      .then(res => res.text())
      .then(text => {
        setContent(text);
        console.log("Loaded markdown:", text);
      });
  }, [isTerms]);

  return (
    <div className="flex justify-center items-stretch min-h-screen w-full" style={{ background: 'linear-gradient(135deg, #101614 0%, #1a2a1a 60%, #1a2320 100%)' }}>
      <div className="w-full max-w-3xl h-full min-h-screen p-8 text-gray-200">
        {content ? (
          <ReactMarkdown
            components={{
              h1: ({children}) => <h1 className="text-2xl font-bold mb-6 text-gray-100">{children}</h1>,
              h2: ({children}) => <h2 className="text-xl font-semibold mt-8 mb-4 text-gray-100">{children}</h2>,
            }}
          >
            {content}
          </ReactMarkdown>
        ) : (
          <div>Loading or no content</div>
        )}
      </div>
    </div>
  );
};

export default TermsPrivacy; 