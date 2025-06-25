import React, { createContext, useContext, useState, useCallback } from "react";

export type ChatSessionContextType = {
  chatSessions: any[];
  selectedSessionId: string | null;
  setSelectedSessionId: (id: string | null) => void;
  refetchChatSessions: () => Promise<void>;
};

const ChatSessionContext = createContext<ChatSessionContextType | undefined>(undefined);

export const ChatSessionProvider: React.FC<{ children: React.ReactNode }> = ({ children }) => {
  console.log("ChatSessionProvider mounted");
  const [chatSessions, setChatSessions] = useState<any[]>([]);
  const [selectedSessionId, setSelectedSessionId] = useState<string | null>(null);

  const refetchChatSessions = useCallback(async () => {
    const res = await fetch(`/api/chat_sessions`, { credentials: 'include' });
    const data = await res.json();
    setChatSessions(data.sessions || []);
    if (data.sessions && data.sessions.length > 0) {
      setSelectedSessionId(data.sessions[0].session_id);
    } else {
      setSelectedSessionId(null);
    }
  }, []);

  return (
    <ChatSessionContext.Provider value={{
      chatSessions,
      selectedSessionId,
      setSelectedSessionId,
      refetchChatSessions
    }}>
      {children}
    </ChatSessionContext.Provider>
  );
};

export const useChatSession = () => {
  const context = useContext(ChatSessionContext);
  if (!context) throw new Error("useChatSession must be used within a ChatSessionProvider");
  return context;
}; 