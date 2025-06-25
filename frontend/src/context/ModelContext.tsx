import React, { createContext, useContext, useState } from "react";

type ModelContextType = {
  selectedModel: string;
  setSelectedModel: (model: string) => void;
};

const ModelContext = createContext<ModelContextType | undefined>(undefined);

export const ModelProvider: React.FC<{ children: React.ReactNode }> = ({ children }) => {
  console.log("ModelProvider mounted");
  const [selectedModel, setSelectedModel] = useState("toxindex-rap");
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