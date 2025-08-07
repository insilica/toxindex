import React, { createContext, useContext, useState, useEffect } from "react";

type ModelContextType = {
  selectedModel: string;
  setSelectedModel: (model: string) => void;
};

const ModelContext = createContext<ModelContextType | undefined>(undefined);

export const ModelProvider: React.FC<{ children: React.ReactNode }> = ({ children }) => {
  console.log('[ModelProvider] ModelProvider mounting');
  const [selectedModel, setSelectedModel] = useState("toxindex-rap");
  
  useEffect(() => {
    console.log('[ModelProvider] ModelProvider mounted');
    return () => {
      console.log('[ModelProvider] ModelProvider unmounting');
    };
  }, []);
  
  console.log('[ModelProvider] ModelProvider rendering with selectedModel:', selectedModel);
  
  return (
    <ModelContext.Provider value={{ selectedModel, setSelectedModel }}>
      {children}
    </ModelContext.Provider>
  );
};

export const useModel = () => {
  const context = useContext(ModelContext);
  if (!context) throw new Error("useModel must be used within a ModelProvider");
  return context;
}; 